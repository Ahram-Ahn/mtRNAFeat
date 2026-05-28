"""Synonymous-recoding thermodynamic permutation test.

Lifted in spirit from `legacy/base_substitution/03.mito_cox1_evo_thermo_pipeline_transcript_span300.py`.

For each (species, gene), takes the wild-type CDS, generates N alternative
sequences from each of three null pools, folds every variant under plain
Vienna MFE (with cfg.max_bp_span as the only constraint), and asks: is the
wild-type ΔG significantly lower (more stable) than each null pool?

**ΔG provenance (this stage NEVER uses the `.db` header MFE value).**
Both reference ΔGs are recomputed from scratch by ViennaRNA on the
sequences/structures in the `.db` records:

  - `WT_MFE`              : ΔG of the wild-type sequence as folded by
                              Vienna MFE (`thermo.fold_mfe`). Apples-to-
                              apples versus the null pools, which are
                              also folded with `thermo.fold_mfe`.
  - `WT_DMS_Eval`         : ΔG of the experimental DMS dot-bracket
                              (read from the `.db` file's structure line)
                              evaluated under Vienna's NN energy model
                              via `thermo.eval_structure`. Asks "is the
                              actual realized in-vivo structure more /
                              less stable than what synonymous shuffling
                              can MFE-fold?"
  - `DMS_FullLen_Vienna`  : same as `WT_DMS_Eval` but for the full
                              CDS DMS structure (no codon truncation),
                              reported once per gene as a sanity check
                              that the chunk-truncated value tracks the
                              full-CDS one. This column is `NaN` when
                              the gene has no DMS structure recorded.

The `.db` header (e.g. `>COX1: -150.4 kcal/mol`) is parsed only by
`io/db_parser.py` and consumed only by reporting stages (`stats`,
`landscape`); it is intentionally bypassed here so substitution ΔGs are
fully reproducible from the dot-bracket structure plus Vienna.

The summary table reports Z-scores and empirical p-values for both refs.

Three nulls (the GC-only ``flat_gc`` and ``positional_gc`` pools were
dropped because the ACGU variants strictly generalize them):
  * **flat_acgu**   — random sequence preserving the gene's full A/C/G/U
                       frequency vector independently. Mitochondrial mRNAs
                       commonly have G≫C on the L-strand (e.g. human ND6:
                       G=191 vs C=37); a symmetric overall-GC null hides
                       that asymmetry entirely.
  * **positional_acgu** — random sequence preserving codon-position
                       A/C/G/U frequencies independently (twelve
                       parameters: 4 nucleotides × 3 codon positions).
                       The right comparator for transcripts where G/C
                       and A/U asymmetries differ across codon positions.
  * **synonymous**   — codon-by-codon synonymous resampling, weighted by
                       observed codon usage AND positional GC. Tests
                       "is structure driven by synonymous-codon choice?"

Default chunk = full CDS up to cfg.substitution_max_nt (300 nt by default —
publishable in a single panel and compatible with mt-mRNA gene lengths).
The 500-nt chunks the user mentioned proved too aggressive for several
short genes (ATP8, ND4L); 300 nt with codon-complete truncation is the
sweet spot. Adjustable via cfg.substitution_max_nt.

Input is CDS-only — sliced via ``io.annotations`` so UTR nucleotides
never enter the codon-aware shuffling (UTRs have no codon table and
would corrupt every codon-aware null).
"""
from __future__ import annotations

import random
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass

import numpy as np
import pandas as pd

from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core import thermo
from mtrnafeat.core.projection import project_structure_to_window, truncate_prefix
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.codons import codon_table_for
from mtrnafeat.io.db_parser import parse_db
from mtrnafeat.progress import progress, step

_SPECIES_SEED_OFFSETS = {
    "human": 101,
    "yeast": 211,
}


@dataclass(frozen=True)
class _PrefixJob:
    species: str
    gene: str
    seq_dna: str
    dms_structure: str
    full_seq_dna: str
    full_dms_structure: str
    n_simulations: int
    max_bp_span: int
    seed: int


def _rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def _dna(seq: str) -> str:
    return seq.upper().replace("U", "T")


def _stable_species_seed_offset(species: str) -> int:
    """Stable per-species seed offset.

    Python's built-in ``hash()`` is intentionally randomized between
    interpreter processes, so it must not be used for reproducible outputs.
    """
    key = species.strip().lower()
    if key in _SPECIES_SEED_OFFSETS:
        return _SPECIES_SEED_OFFSETS[key]
    return sum((i + 1) * ord(ch) for i, ch in enumerate(key)) % 997


def _truncate_codon_aligned(seq: str, max_nt: int | None) -> str:
    seq = _dna(seq)
    if max_nt is None or len(seq) <= max_nt:
        n = len(seq) - (len(seq) % 3)
        return seq[:n]
    return seq[: (max_nt - max_nt % 3)]


def _split_codons(seq_dna: str) -> list[str]:
    return [seq_dna[i:i + 3] for i in range(0, len(seq_dna), 3)]


def _build_synonyms(table: dict[str, str]) -> dict[str, list[str]]:
    syns: dict[str, list[str]] = defaultdict(list)
    for codon, aa in table.items():
        syns[aa].append(codon)
    return dict(syns)


def _gc_fraction(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / max(1, len(s))


def _acgu_fraction(seq_dna: str) -> dict[str, float]:
    """Return the gene's overall A / C / G / T(=U) frequency vector.

    Frequencies are normalized over A+C+G+T only (any non-ACGT character
    is ignored, which makes the function robust to ambiguity codes).
    """
    s = seq_dna.upper()
    counts = {nt: s.count(nt) for nt in "ACGT"}
    total = sum(counts.values()) or 1
    return {nt: counts[nt] / total for nt in "ACGT"}


def _positional_gc(seq_dna: str) -> dict[int, float]:
    codons = _split_codons(seq_dna)
    counts = {1: [0, 0], 2: [0, 0], 3: [0, 0]}  # [gc, total]
    for c in codons:
        if len(c) != 3:
            continue
        for pos, nt in enumerate(c, start=1):
            counts[pos][1] += 1
            if nt in "GC":
                counts[pos][0] += 1
    return {p: (gc / tot) if tot else 0.0 for p, (gc, tot) in counts.items()}


def _positional_acgu(seq_dna: str) -> dict[int, dict[str, float]]:
    """Per-codon-position A / C / G / T frequency vectors."""
    codons = _split_codons(seq_dna)
    counts: dict[int, dict[str, int]] = {p: {nt: 0 for nt in "ACGT"} for p in (1, 2, 3)}
    totals: dict[int, int] = {1: 0, 2: 0, 3: 0}
    for c in codons:
        if len(c) != 3:
            continue
        for pos, nt in enumerate(c, start=1):
            if nt in counts[pos]:
                counts[pos][nt] += 1
                totals[pos] += 1
    out: dict[int, dict[str, float]] = {}
    for pos in (1, 2, 3):
        denom = totals[pos] or 1
        out[pos] = {nt: counts[pos][nt] / denom for nt in "ACGT"}
    return out


def _codon_usage_prior(seq_dna: str, table: dict[str, str]) -> dict[str, dict[str, float]]:
    syns = _build_synonyms(table)
    counts: dict[str, Counter] = defaultdict(Counter)
    for c in _split_codons(seq_dna):
        if c in table:
            counts[table[c]][c] += 1
    prior: dict[str, dict[str, float]] = {}
    for aa, options in syns.items():
        n_obs = counts[aa]
        total = sum(n_obs.values()) + len(options)  # Laplace smoothing
        prior[aa] = {opt: (n_obs[opt] + 1) / total for opt in options}
    return prior


# --- pool samplers ---------------------------------------------------------

def _flat_gc(length: int, gc: float, rng: random.Random) -> str:
    out = []
    for _ in range(length):
        if rng.random() < gc:
            out.append(rng.choice("GC"))
        else:
            out.append(rng.choice("AT"))
    return "".join(out)


def _flat_acgu(length: int, freqs: dict[str, float], rng: random.Random) -> str:
    """IID nucleotide draws preserving the full A/C/G/T frequency vector.

    Uses ``random.choices`` with explicit weights so G≠C and A≠T
    asymmetries are respected. Reduces to ``_flat_gc`` only when
    P(G)=P(C) and P(A)=P(T).
    """
    nts = "ACGT"
    weights = [freqs.get(nt, 0.0) for nt in nts]
    return "".join(rng.choices(nts, weights=weights, k=length))


def _positional_gc_sample(n_codons: int, pos_gc: dict[int, float], rng: random.Random) -> str:
    out = []
    for _ in range(n_codons):
        for pos in (1, 2, 3):
            p = pos_gc[pos]
            out.append(rng.choice("GC") if rng.random() < p else rng.choice("AT"))
    return "".join(out)


def _positional_acgu_sample(n_codons: int,
                            pos_freqs: dict[int, dict[str, float]],
                            rng: random.Random) -> str:
    """Per-codon-position 4-way A/C/G/T sampling.

    Generalizes ``_positional_gc_sample`` by tracking each nucleotide's
    frequency at each codon position independently (12 parameters total).
    """
    nts = "ACGT"
    out: list[str] = []
    for _ in range(n_codons):
        for pos in (1, 2, 3):
            freqs = pos_freqs[pos]
            weights = [freqs.get(nt, 0.0) for nt in nts]
            out.append(rng.choices(nts, weights=weights, k=1)[0])
    return "".join(out)


def _synonymous_sample(seq_dna: str, table: dict[str, str], pos_gc: dict[int, float],
                        usage_prior: dict[str, dict[str, float]], rng: random.Random) -> str:
    syns = _build_synonyms(table)
    out: list[str] = []
    for c in _split_codons(seq_dna):
        if c not in table or table[c] in {"*"}:
            out.append(c)
            continue
        aa = table[c]
        candidates = syns[aa]
        weights = []
        for cand in candidates:
            w = usage_prior.get(aa, {}).get(cand, 1.0 / len(candidates))
            for pos, nt in enumerate(cand, start=1):
                p = pos_gc[pos]
                w *= p if nt in "GC" else (1.0 - p)
            weights.append(max(1e-9, w))
        total = sum(weights)
        threshold = rng.random() * total
        cumulative = 0.0
        chosen = candidates[-1]
        for cand, w in zip(candidates, weights, strict=True):
            cumulative += w
            if cumulative >= threshold:
                chosen = cand
                break
        out.append(chosen)
    return "".join(out)


# --- per-gene runner -------------------------------------------------------

def _run_one(job: _PrefixJob) -> pd.DataFrame:
    rng = random.Random(job.seed)
    seq_dna = _truncate_codon_aligned(job.seq_dna, max_nt=None)
    if not seq_dna:
        return pd.DataFrame()
    table = codon_table_for(job.species)
    pos_gc = _positional_gc(seq_dna)
    pos_acgu = _positional_acgu(seq_dna)
    usage_prior = _codon_usage_prior(seq_dna, table)
    overall_acgu = _acgu_fraction(seq_dna)
    n_codons = len(seq_dna) // 3
    length = len(seq_dna)

    # Wild-type ΔG references — both recomputed by ViennaRNA on the .db's
    # sequence/structure (the .db header MFE is intentionally NOT used):
    #   wt_mfe           : Vienna MFE of WT sequence (apples-to-apples vs pool)
    #   wt_dms           : Vienna eval of the .db dot-bracket TRUNCATED to
    #                      the chunk length being permuted
    #   wt_dms_fulllen   : Vienna eval of the .db dot-bracket at the FULL
    #                      transcript length (sanity check vs. wt_dms; not
    #                      compared to the pool)
    rna_wt = _rna(seq_dna)
    _, wt_mfe = thermo.fold_mfe(rna_wt, max_bp_span=job.max_bp_span)
    wt_dms = float("nan")
    if job.dms_structure:
        dms_pref = truncate_prefix(job.dms_structure, length)
        if dms_pref:
            try:
                wt_dms = thermo.eval_structure(rna_wt, dms_pref, max_bp_span=job.max_bp_span)
            except Exception:
                wt_dms = float("nan")
    wt_dms_fulllen = float("nan")
    if job.full_dms_structure and job.full_seq_dna:
        try:
            wt_dms_fulllen = thermo.eval_structure(
                _rna(job.full_seq_dna), job.full_dms_structure,
                max_bp_span=job.max_bp_span,
            )
        except Exception:
            wt_dms_fulllen = float("nan")

    # GC-only pools (flat_gc, positional_gc) were dropped at user request —
    # the ACGU-aware versions strictly generalize them (they reduce to the
    # GC variants when G=C and A=U) and the symmetric-GC nulls hide the
    # H-strand C/G and A/U asymmetries we care about.
    pools = {
        "synonymous": [_synonymous_sample(seq_dna, table, pos_gc, usage_prior, rng)
                       for _ in range(job.n_simulations)],
        "flat_acgu": [_flat_acgu(length, overall_acgu, rng)
                      for _ in range(job.n_simulations)],
        "positional_acgu": [_positional_acgu_sample(n_codons, pos_acgu, rng)
                            for _ in range(job.n_simulations)],
    }
    rows = [
        {
            "Species": job.species,
            "Gene": canonical_gene(job.gene),
            "Pool": "WildType_MFE",
            "Simulation": 0,
            "MFE_kcal_per_mol": float(wt_mfe),
            "MFE_kcal_per_nt": float(wt_mfe) / length if length else float("nan"),
            "Length_nt": length,
            "Sequence_Scope": "CDS codon-aligned prefix",
            "Max_BP_Span": int(job.max_bp_span),
        },
        {
            "Species": job.species,
            "Gene": canonical_gene(job.gene),
            "Pool": "WildType_DMS_Eval",
            "Simulation": 0,
            "MFE_kcal_per_mol": float(wt_dms),
            "MFE_kcal_per_nt": float(wt_dms) / length if length else float("nan"),
            "Length_nt": length,
            "Sequence_Scope": "CDS codon-aligned prefix",
            "Max_BP_Span": int(job.max_bp_span),
        },
        {
            "Species": job.species,
            "Gene": canonical_gene(job.gene),
            "Pool": "WildType_DMS_Eval_FullLength",
            "Simulation": 0,
            "MFE_kcal_per_mol": float(wt_dms_fulllen),
            "MFE_kcal_per_nt": (
                float(wt_dms_fulllen) / len(job.full_seq_dna)
                if job.full_seq_dna else float("nan")
            ),
            "Length_nt": len(job.full_seq_dna),
            "Sequence_Scope": "full CDS",
            "Max_BP_Span": int(job.max_bp_span),
        },
    ]
    for pool_name, seqs in pools.items():
        for i, s in enumerate(seqs, start=1):
            _, dg = thermo.fold_mfe(_rna(s), max_bp_span=job.max_bp_span)
            rows.append({
                "Species": job.species,
                "Gene": canonical_gene(job.gene),
                "Pool": pool_name,
                "Simulation": i,
                "MFE_kcal_per_mol": float(dg),
                "MFE_kcal_per_nt": float(dg) / length if length else float("nan"),
                "Length_nt": length,
                "Sequence_Scope": "CDS codon-aligned prefix",
                "Max_BP_Span": int(job.max_bp_span),
            })
    return pd.DataFrame(rows)


def _summarize(dist: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (species, gene), g in dist.groupby(["Species", "Gene"]):
        wt_mfe_row = g[g["Pool"] == "WildType_MFE"]["MFE_kcal_per_mol"]
        wt_dms_row = g[g["Pool"] == "WildType_DMS_Eval"]["MFE_kcal_per_mol"]
        wt_dms_full_row = g[g["Pool"] == "WildType_DMS_Eval_FullLength"]["MFE_kcal_per_mol"]
        if wt_mfe_row.empty:
            continue
        wt_mfe = float(wt_mfe_row.iloc[0])
        wt_dms = float(wt_dms_row.iloc[0]) if not wt_dms_row.empty else float("nan")
        wt_dms_full = float(wt_dms_full_row.iloc[0]) if not wt_dms_full_row.empty else float("nan")
        chunk_len = int(g[g["Pool"] == "WildType_MFE"]["Length_nt"].iloc[0])
        full_len_row = g[g["Pool"] == "WildType_DMS_Eval_FullLength"]["Length_nt"]
        full_len = int(full_len_row.iloc[0]) if not full_len_row.empty else 0
        for pool in ("flat_acgu", "positional_acgu", "synonymous"):
            pool_vals = g[g["Pool"] == pool]["MFE_kcal_per_mol"].values.astype(float)
            if len(pool_vals) == 0:
                continue
            mean = float(np.mean(pool_vals))
            sd = float(np.std(pool_vals, ddof=1)) if len(pool_vals) > 1 else 0.0
            n = len(pool_vals)

            def _z_p(
                observed: float,
                vals: np.ndarray = pool_vals,
                pool_mean: float = mean,
                pool_sd: float = sd,
                pool_n: int = n,
            ) -> tuple[float, float]:
                if not np.isfinite(observed):
                    return float("nan"), float("nan")
                z = (observed - pool_mean) / pool_sd if pool_sd > 0 else 0.0
                p_lower = (1 + int(np.sum(vals <= observed))) / (pool_n + 1)
                return float(z), float(p_lower)

            z_mfe, p_mfe = _z_p(wt_mfe)
            z_dms, p_dms = _z_p(wt_dms)
            delta_mfe = wt_mfe - mean
            delta_dms = wt_dms - mean if np.isfinite(wt_dms) else float("nan")
            rows.append({
                "Species": species,
                "Gene": gene,
                "Pool": pool,
                "Length_nt": chunk_len,
                "DMS_Full_CDS_Length_nt": full_len,
                "WT_MFE": wt_mfe,
                "WT_MFE_per_nt": wt_mfe / chunk_len if chunk_len else float("nan"),
                "WT_DMS_Eval": wt_dms,
                "WT_DMS_Eval_per_nt": wt_dms / chunk_len if chunk_len else float("nan"),
                "WT_DMS_Eval_FullLength": wt_dms_full,
                "WT_DMS_Eval_FullLength_per_nt": (
                    wt_dms_full / full_len if full_len else float("nan")
                ),
                "Pool_Mean_MFE": mean,
                "Pool_Mean_MFE_per_nt": mean / chunk_len if chunk_len else float("nan"),
                "Pool_SD_MFE": sd,
                "Delta_WT_MFE_minus_Pool_Mean": delta_mfe,
                "Delta_WT_MFE_per_nt_minus_Pool_Mean": (
                    delta_mfe / chunk_len if chunk_len else float("nan")
                ),
                "Delta_WT_DMS_minus_Pool_Mean": delta_dms,
                "Delta_WT_DMS_per_nt_minus_Pool_Mean": (
                    delta_dms / chunk_len if chunk_len and np.isfinite(delta_dms) else float("nan")
                ),
                "Z_WT_MFE_vs_Pool": z_mfe,
                "Empirical_p_WT_MFE_more_stable": p_mfe,
                "Pool_Percentile_WT_MFE": 100.0 * p_mfe if np.isfinite(p_mfe) else float("nan"),
                "Z_WT_DMS_vs_Pool": z_dms,
                "Empirical_p_WT_DMS_more_stable": p_dms,
                "Pool_Percentile_WT_DMS": 100.0 * p_dms if np.isfinite(p_dms) else float("nan"),
                "N_Simulations": n,
                "DMS_dG_Source": (
                    "Vienna eval_structure on projected .db CDS dot-bracket; "
                    "WT_DMS_Eval is the codon-aligned CDS prefix, "
                    "WT_DMS_Eval_FullLength is the full CDS"
                ),
            })
    return pd.DataFrame(rows)


def _collect_sequential(jobs: list[_PrefixJob]) -> list[pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    for j in progress(jobs, desc="substitution (genes)", unit="gene"):
        frames.append(_run_one(j))
    return frames


def run_substitution_thermo(cfg: Config) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (long_distribution_df, summary_df)."""
    step(f"running substitution-thermo (n={cfg.substitution_n_simulations} per pool, Vienna MFE)")
    jobs: list[_PrefixJob] = []
    rng_seed_base = int(cfg.seed)
    for species, fname in cfg.db_files.items():
        rec_by_gene = {r.gene: r for r in parse_db(cfg.data_dir / fname)}
        for gi, gene in enumerate(cfg.target_genes):
            target = canonical_gene(gene)
            if target not in rec_by_gene:
                continue
            rec = rec_by_gene[target]

            # Synonymous-codon shuffling requires CDS-only input — UTR
            # nucleotides have no codon table and would corrupt every
            # codon-aware null. Slice using the species annotation; the
            # structure must be projected (not just sliced) so any pairs
            # crossing the CDS boundary become '.' instead of leaving
            # an unbalanced dot-bracket that segfaults ViennaRNA.
            try:
                annot = annotation_for(species, target)
                l_utr5 = int(annot["l_utr5"])
                l_cds = int(annot["l_cds"])
                cds_seq = rec.sequence[l_utr5:l_utr5 + l_cds]
                if rec.structure:
                    cds_struct = project_structure_to_window(
                        rec.structure, l_utr5, l_utr5 + l_cds,
                    )
                else:
                    cds_struct = ""
            except KeyError:
                cds_seq = rec.sequence
                cds_struct = rec.structure or ""

            seq = _truncate_codon_aligned(cds_seq, max_nt=cfg.substitution_max_nt)
            if len(seq) < 60:
                continue
            full_seq = _dna(cds_seq)
            full_struct = cds_struct
            jobs.append(_PrefixJob(
                species=species, gene=target, seq_dna=seq,
                dms_structure=cds_struct,
                full_seq_dna=full_seq,
                full_dms_structure=full_struct,
                n_simulations=int(cfg.substitution_n_simulations),
                max_bp_span=int(cfg.max_bp_span),
                seed=rng_seed_base + 1000 * gi + _stable_species_seed_offset(species),
            ))

    if not jobs:
        return pd.DataFrame(), pd.DataFrame()

    workers = max(1, int(cfg.n_workers))
    if workers == 1 or len(jobs) == 1:
        frames = _collect_sequential(jobs)
    else:
        try:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                futures = {ex.submit(_run_one, j): j for j in jobs}
                frames = []
                for fut in progress(as_completed(futures), desc="substitution (genes)",
                                      total=len(futures), unit="gene"):
                    frames.append(fut.result())
        except (OSError, PermissionError) as exc:
            step(
                "substitution multiprocessing unavailable "
                f"({type(exc).__name__}: {exc}); falling back to sequential execution"
            )
            frames = _collect_sequential(jobs)

    dist = pd.concat([f for f in frames if not f.empty], ignore_index=True) if frames else pd.DataFrame()
    summary = _summarize(dist) if not dist.empty else pd.DataFrame()
    return dist, summary
