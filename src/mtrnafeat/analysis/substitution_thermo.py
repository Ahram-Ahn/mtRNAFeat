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
  - `DMS_FullLen_Vienna`  : same as `WT_DMS_Eval` but for the FULL
                              transcript-length DMS structure (no codon
                              truncation), reported once per gene as a
                              sanity check that the chunk-truncated
                              value tracks the full-length one. This
                              column is `NaN` when the gene has no DMS
                              structure recorded.

The `.db` header (e.g. `>COX1: -150.4 kcal/mol`) is parsed only by
`io/db_parser.py` and consumed only by reporting stages (`stats`,
`landscape`); it is intentionally bypassed here so substitution ΔGs are
fully reproducible from the dot-bracket structure plus Vienna.

The summary table reports Z-scores and empirical p-values for both refs.

Three nulls:
  * **flat_gc**     — random sequence at the gene's overall GC content;
                       AT and GC drawn IID. Tests "is structure driven by GC?"
  * **positional_gc** — random sequence preserving codon-position GC content
                       (1st, 2nd, wobble each fixed). Tests "is structure driven
                       by codon-bias of the genome?"
  * **synonymous**   — codon-by-codon synonymous resampling, weighted by
                       observed codon usage AND positional GC. Tests
                       "is structure driven by synonymous-codon choice?"

Default chunk = full CDS up to cfg.substitution_max_nt (300 nt by default —
publishable in a single panel and compatible with mt-mRNA gene lengths).
The 500-nt chunks the user mentioned proved too aggressive for several
short genes (ATP8, ND4L); 300 nt with codon-complete truncation is the
sweet spot. Adjustable via cfg.substitution_max_nt.
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
from mtrnafeat.core.projection import truncate_prefix
from mtrnafeat.io.codons import codon_table_for
from mtrnafeat.io.db_parser import parse_db
from mtrnafeat.progress import progress, step


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


def _positional_gc_sample(n_codons: int, pos_gc: dict[int, float], rng: random.Random) -> str:
    out = []
    for _ in range(n_codons):
        for pos in (1, 2, 3):
            p = pos_gc[pos]
            out.append(rng.choice("GC") if rng.random() < p else rng.choice("AT"))
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
        for cand, w in zip(candidates, weights):
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
    usage_prior = _codon_usage_prior(seq_dna, table)
    overall_gc = _gc_fraction(seq_dna)
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

    pools = {
        "flat_gc": [_flat_gc(length, overall_gc, rng) for _ in range(job.n_simulations)],
        "positional_gc": [_positional_gc_sample(n_codons, pos_gc, rng)
                          for _ in range(job.n_simulations)],
        "synonymous": [_synonymous_sample(seq_dna, table, pos_gc, usage_prior, rng)
                       for _ in range(job.n_simulations)],
    }
    rows = [
        {
            "Species": job.species,
            "Gene": canonical_gene(job.gene),
            "Pool": "WildType_MFE",
            "Simulation": 0,
            "MFE_kcal_per_mol": float(wt_mfe),
            "Length_nt": length,
        },
        {
            "Species": job.species,
            "Gene": canonical_gene(job.gene),
            "Pool": "WildType_DMS_Eval",
            "Simulation": 0,
            "MFE_kcal_per_mol": float(wt_dms),
            "Length_nt": length,
        },
        {
            "Species": job.species,
            "Gene": canonical_gene(job.gene),
            "Pool": "WildType_DMS_Eval_FullLength",
            "Simulation": 0,
            "MFE_kcal_per_mol": float(wt_dms_fulllen),
            "Length_nt": len(job.full_seq_dna),
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
                "Length_nt": length,
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
        for pool in ("flat_gc", "positional_gc", "synonymous"):
            pool_vals = g[g["Pool"] == pool]["MFE_kcal_per_mol"].values.astype(float)
            if len(pool_vals) == 0:
                continue
            mean = float(np.mean(pool_vals))
            sd = float(np.std(pool_vals, ddof=1)) if len(pool_vals) > 1 else 0.0
            n = len(pool_vals)

            def _z_p(observed: float) -> tuple[float, float]:
                if not np.isfinite(observed):
                    return float("nan"), float("nan")
                z = (observed - mean) / sd if sd > 0 else 0.0
                p_lower = (1 + int(np.sum(pool_vals <= observed))) / (n + 1)
                return float(z), float(p_lower)

            z_mfe, p_mfe = _z_p(wt_mfe)
            z_dms, p_dms = _z_p(wt_dms)
            rows.append({
                "Species": species,
                "Gene": gene,
                "Pool": pool,
                "Length_nt": int(g[g["Pool"] == "WildType_MFE"]["Length_nt"].iloc[0]),
                "WT_MFE": wt_mfe,
                "WT_DMS_Eval": wt_dms,
                "WT_DMS_Eval_FullLength": wt_dms_full,
                "Pool_Mean_MFE": mean,
                "Pool_SD_MFE": sd,
                "Z_WT_MFE_vs_Pool": z_mfe,
                "Empirical_p_WT_MFE_more_stable": p_mfe,
                "Z_WT_DMS_vs_Pool": z_dms,
                "Empirical_p_WT_DMS_more_stable": p_dms,
                "N_Simulations": n,
                "DMS_dG_Source": "Vienna eval_structure on .db dot-bracket",
            })
    return pd.DataFrame(rows)


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
            seq = _truncate_codon_aligned(rec.sequence, max_nt=cfg.substitution_max_nt)
            if len(seq) < 60:
                continue
            full_seq = _dna(rec.sequence)
            full_struct = rec.structure or ""
            jobs.append(_PrefixJob(
                species=species, gene=target, seq_dna=seq,
                dms_structure=rec.structure or "",
                full_seq_dna=full_seq,
                full_dms_structure=full_struct,
                n_simulations=int(cfg.substitution_n_simulations),
                max_bp_span=int(cfg.max_bp_span),
                seed=rng_seed_base + 1000 * gi + hash(species) % 997,
            ))

    if not jobs:
        return pd.DataFrame(), pd.DataFrame()

    workers = max(1, int(cfg.n_workers))
    frames: list[pd.DataFrame] = []
    if workers == 1 or len(jobs) == 1:
        for j in progress(jobs, desc="substitution (genes)", unit="gene"):
            frames.append(_run_one(j))
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(_run_one, j): j for j in jobs}
            for fut in progress(as_completed(futures), desc="substitution (genes)",
                                  total=len(futures), unit="gene"):
                frames.append(fut.result())

    dist = pd.concat([f for f in frames if not f.empty], ignore_index=True) if frames else pd.DataFrame()
    summary = _summarize(dist) if not dist.empty else pd.DataFrame()
    return dist, summary
