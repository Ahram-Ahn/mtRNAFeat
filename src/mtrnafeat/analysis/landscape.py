"""GC-gradient simulation + experimental overlay.

Replaces legacy 01, 02, 03.test_02, and 03_integrate_sim_dms.

Three products:
- `simulate_specific_conditions`: uses **empirical per-species (A, U, G, C)
  frequencies** by default (see `compute_empirical_freqs`). Falls back to
  `Config.sim_gc_conditions` (symmetric GC) only if a species' empirical
  frequencies are unavailable. This is the change requested by the user:
  human mt-mRNAs are C-enriched on H-strand transcripts (C ≫ G), and a
  symmetric GC null misses that.
- `simulate_gradient`: continuous 0%–100% GC sweep (still symmetric — by
  design, since the question of the gradient is "what does GC alone do?").
- `experimental_overlay`: per-transcript (foldedness, normMFE, paired-pair-type)
  loaded from .db files.

Progress bars via mtrnafeat.progress.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from mtrnafeat.analysis.statistics import paired_composition, sequence_gc_pct
from mtrnafeat.config import Config
from mtrnafeat.core import thermo
from mtrnafeat.core.projection import project_structure_to_window
from mtrnafeat.core.shuffle import random_gc_sequence, random_sequence_with_freqs
from mtrnafeat.core.structure import paired_fraction
from mtrnafeat.io.annotations import annotation_for_sample, infer_annotation_species
from mtrnafeat.io.db_parser import parse_db
from mtrnafeat.progress import progress, step
from mtrnafeat.rng import make_rng

#: Genes to drop when computing **heavy-strand-only** species frequencies.
#: Human ND6 is the sole H-strand-template / L-strand-encoded mt-mRNA, so its
#: base composition is the opposite skew (G≫C, U≫A) of the other 12 H-strand
#: mRNAs and would dilute the H-strand C-enrichment if included.
HEAVY_STRAND_EXCLUSIONS: dict[str, frozenset[str]] = {
    "Human": frozenset({"ND6"}),
}


def compute_empirical_freqs(db_path, exclude_genes: frozenset[str] = frozenset()) -> dict[str, float]:
    """Aggregate (A, U, G, C) frequency over records in a .db file.

    If ``exclude_genes`` is non-empty, records whose canonical gene name is in
    that set are skipped — used to compute heavy-strand-only frequencies by
    dropping L-strand-encoded transcripts (e.g. human ND6).
    """
    counts: Counter[str] = Counter()
    for rec in parse_db(db_path):
        if rec.gene in exclude_genes:
            continue
        counts.update(rec.sequence.upper())
    total = sum(counts[b] for b in "AUGC")
    if total == 0:
        return {"A": 0.25, "U": 0.25, "G": 0.25, "C": 0.25}
    return {b: counts[b] / total for b in "AUGC"}


def _base_composition(seq: str) -> dict[str, float]:
    s = seq.upper()
    n = len(s)
    if n == 0:
        return {b: 0.0 for b in "ACGU"}
    return {b: 100.0 * s.count(b) / n for b in "ACGU"}


def species_freqs_for_pipeline(cfg: Config, heavy_strand: bool = False) -> dict[str, dict[str, float]]:
    """Return {species: {A,U,G,C}} for every species declared in cfg.db_files,
    using either cfg.sim_freqs_per_species (if provided) or empirical .db data.

    When ``heavy_strand`` is True, empirical frequencies are computed from
    H-strand-encoded transcripts only (see :data:`HEAVY_STRAND_EXCLUSIONS`).
    Species with no exclusions defined behave identically to the default.
    Explicit overrides in ``cfg.sim_freqs_per_species`` still win.
    """
    out: dict[str, dict[str, float]] = {}
    given = getattr(cfg, "sim_freqs_per_species", None) or {}
    for species, fname in cfg.db_files.items():
        if species in given:
            out[species] = {b: float(given[species].get(b, 0.0)) for b in "AUGC"}
        else:
            exclude = HEAVY_STRAND_EXCLUSIONS.get(species, frozenset()) if heavy_strand else frozenset()
            out[species] = compute_empirical_freqs(cfg.data_dir / fname, exclude_genes=exclude)
    return out


def simulate_specific_conditions(cfg: Config) -> pd.DataFrame:
    """Simulate `Sim {Species}` clouds using per-species empirical freqs.

    Adds two extra clouds for the published symmetric-GC nulls (Sim Yeast 5'UTR
    7%, Sim Yeast CDS 30%) so the legacy comparison panels still render.
    """
    rng = make_rng(cfg.seed)
    rows: list[dict] = []
    species_freqs = species_freqs_for_pipeline(cfg)
    step(f"simulating empirical-frequency clouds for: {', '.join(species_freqs)}")

    # Empirical per-species clouds.
    for species, freqs in species_freqs.items():
        label = (f"Sim {species} ({freqs['A']:.2f}A {freqs['U']:.2f}U "
                 f"{freqs['G']:.2f}G {freqs['C']:.2f}C)")
        for _ in progress(range(cfg.sim_num_sequences),
                           desc=f"sim {species} (empirical)", unit="seq"):
            seq = random_sequence_with_freqs(cfg.sim_seq_length, freqs, rng)
            struct, mfe = thermo.fold_mfe(seq)
            length = len(seq)
            gc_pair, au_pair, gu_pair = paired_composition(seq, struct)
            rows.append({
                "Condition": label, "Data_Type": "Simulation",
                "Gene": "Simulated", "Species": species, "Length": length,
                "MFE": mfe, "Normalized_MFE_per_nt": mfe / length,
                "Foldedness_Pct": 100.0 * paired_fraction(struct),
                "Sequence_GC_Pct": sequence_gc_pct(seq),
                "Paired_GC_Pct": gc_pair,
                "Paired_AU_Pct": au_pair,
                "Paired_GU_Pct": gu_pair,
            })

    # Symmetric-GC clouds (legacy reference points, still useful).
    for label, gc in cfg.sim_gc_conditions.items():
        for _ in progress(range(cfg.sim_num_sequences),
                           desc=f"sim {label}", unit="seq"):
            seq = random_gc_sequence(cfg.sim_seq_length, gc, rng)
            struct, mfe = thermo.fold_mfe(seq)
            length = len(seq)
            gc_pair, au_pair, gu_pair = paired_composition(seq, struct)
            rows.append({
                "Condition": label, "Data_Type": "Simulation",
                "Gene": "Simulated", "Species": "n/a", "Length": length,
                "MFE": mfe, "Normalized_MFE_per_nt": mfe / length,
                "Foldedness_Pct": 100.0 * paired_fraction(struct),
                "Sequence_GC_Pct": sequence_gc_pct(seq),
                "Paired_GC_Pct": gc_pair,
                "Paired_AU_Pct": au_pair,
                "Paired_GU_Pct": gu_pair,
            })

    return pd.DataFrame(rows)


def simulate_gradient(cfg: Config) -> pd.DataFrame:
    rng = make_rng(cfg.seed + 1)
    rows: list[dict] = []
    gc_steps = np.linspace(0.0, 1.0, cfg.gradient_steps)
    step(f"simulating GC gradient ({cfg.gradient_steps} steps × {cfg.gradient_seqs_per_step} sequences)")
    for gc in progress(gc_steps, desc="GC gradient", unit="GC"):
        for _ in range(cfg.gradient_seqs_per_step):
            seq = random_gc_sequence(cfg.sim_seq_length, float(gc), rng)
            struct, mfe = thermo.fold_mfe(seq)
            length = len(seq)
            gc_pair, au_pair, gu_pair = paired_composition(seq, struct)
            rows.append({
                "Condition": "Gradient", "Data_Type": "Simulation",
                "Sequence_GC_Pct": sequence_gc_pct(seq),
                "GC_Target_Pct": float(gc) * 100.0,
                "Foldedness_Pct": 100.0 * paired_fraction(struct),
                "Normalized_MFE_per_nt": mfe / length,
                "Paired_GC_Pct": gc_pair,
                "Paired_AU_Pct": au_pair,
                "Paired_GU_Pct": gu_pair,
            })
    return pd.DataFrame(rows)


def experimental_overlay(cfg: Config) -> pd.DataFrame:
    rows: list[dict] = []
    for species, fname in cfg.db_files.items():
        path = cfg.data_dir / fname
        for rec in parse_db(path):
            length = len(rec.sequence)
            try:
                dms_eval = thermo.eval_structure(rec.sequence, rec.structure)
            except Exception:
                dms_eval = float("nan")
            gc_pair, au_pair, gu_pair = paired_composition(rec.sequence, rec.structure)
            base = _base_composition(rec.sequence)
            rows.append({
                "Condition": f"Exp: {species} DMS",
                "Data_Type": "Experimental",
                "Gene": rec.gene, "Species": species, "Length": length,
                "MFE": dms_eval,
                "Header_MFE": rec.mfe,
                "MFE_Source": "Vienna eval_structure(sequence, .db dot-bracket)",
                "Normalized_MFE_per_nt": dms_eval / length if length else 0.0,
                "Foldedness_Pct": 100.0 * paired_fraction(rec.structure),
                "Sequence_GC_Pct": sequence_gc_pct(rec.sequence),
                "Paired_GC_Pct": gc_pair,
                "Paired_AU_Pct": au_pair,
                "Paired_GU_Pct": gu_pair,
                "Pct_A": base["A"],
                "Pct_C": base["C"],
                "Pct_G": base["G"],
                "Pct_U": base["U"],
            })
    return pd.DataFrame(rows)


def experimental_overlay_regions(cfg: Config) -> pd.DataFrame:
    """Region-level DMS points with de novo Vienna ΔG evaluation.

    Rows are emitted only when a sample label can be mapped to a bundled
    annotation species. This currently supports Human/Yeast labels directly
    and arbitrary sample names via ``cfg.sample_annotation_species`` or the
    Human/Yeast token inference in ``io.annotations``.
    """
    rows: list[dict] = []
    mapping = getattr(cfg, "sample_annotation_species", None) or {}
    for species, fname in cfg.db_files.items():
        annotation_species = infer_annotation_species(species, mapping)
        if annotation_species is None:
            continue
        path = cfg.data_dir / fname
        for rec in parse_db(path):
            try:
                annot = annotation_for_sample(species, rec.gene, mapping)
            except KeyError:
                continue
            regions = (
                ("5'UTR", 0, int(annot["l_utr5"])),
                ("CDS", int(annot["l_utr5"]), int(annot["l_utr5"]) + int(annot["l_cds"])),
                (
                    "3'UTR/tail",
                    int(annot["l_utr5"]) + int(annot["l_cds"]),
                    min(len(rec.sequence), int(annot["l_tr"])),
                ),
            )
            for region, start, end in regions:
                start = max(0, min(start, len(rec.sequence)))
                end = max(start, min(end, len(rec.sequence)))
                if end <= start:
                    continue
                seq_w = rec.sequence[start:end]
                struct_w = project_structure_to_window(rec.structure, start, end)
                try:
                    dms_eval = thermo.eval_structure(seq_w, struct_w)
                except Exception:
                    dms_eval = float("nan")
                gc_pair, au_pair, gu_pair = paired_composition(seq_w, struct_w)
                base = _base_composition(seq_w)
                length = len(seq_w)
                rows.append({
                    "Condition": f"Exp: {species} DMS {region}",
                    "Data_Type": "Experimental",
                    "Gene": rec.gene,
                    "Species": species,
                    "Annotation_Species": annotation_species,
                    "Region": region,
                    "Region_Start_1based": start + 1,
                    "Region_End_1based": end,
                    "Length": length,
                    "MFE": dms_eval,
                    "MFE_Source": "Vienna eval_structure(region sequence, projected .db dot-bracket)",
                    "Normalized_MFE_per_nt": dms_eval / length if length else 0.0,
                    "Foldedness_Pct": 100.0 * paired_fraction(struct_w),
                    "Sequence_GC_Pct": sequence_gc_pct(seq_w),
                    "Paired_GC_Pct": gc_pair,
                    "Paired_AU_Pct": au_pair,
                    "Paired_GU_Pct": gu_pair,
                    "Pct_A": base["A"],
                    "Pct_C": base["C"],
                    "Pct_G": base["G"],
                    "Pct_U": base["U"],
                })
    return pd.DataFrame(rows)


def simulate_biased_gradient(cfg: Config, heavy_strand: bool = False) -> pd.DataFrame:
    """Per-species GC gradient that PRESERVES the empirical C:G and A:U ratios.

    The standard ``simulate_gradient`` uses symmetric base ratios (G=C, A=U).
    This variant fixes the species-specific C/(G+C) and A/(A+U) splits across
    the whole 0–100% GC sweep, so the resulting null cloud reflects how
    species-biased strand chemistry (e.g. human heavy strand: C ≫ G, A ≫ U)
    behaves along the same GC axis.

    With ``heavy_strand=True``, species ratios come from H-strand transcripts
    only (drops human ND6, which is L-strand-encoded and reverses the skew).
    The "Condition" label is suffixed with " — H-strand" so downstream code
    can distinguish the two clouds when both are merged.
    """
    seed_offset = 17 if heavy_strand else 7
    rng = make_rng(cfg.seed + seed_offset)
    rows: list[dict] = []
    species_freqs = species_freqs_for_pipeline(cfg, heavy_strand=heavy_strand)
    gc_steps = np.linspace(0.0, 1.0, cfg.gradient_steps)
    label_tag = "biased H-strand GC gradient" if heavy_strand else "biased GC gradient"
    step(f"simulating {label_tag} per species "
         f"({cfg.gradient_steps} steps × {cfg.gradient_seqs_per_step} seq)")

    for species, freqs in species_freqs.items():
        emp_gc = freqs["G"] + freqs["C"]
        emp_au = freqs["A"] + freqs["U"]
        # Split ratios within GC and AU; default to symmetric if a class is
        # empty in the empirical data (unlikely but safe).
        c_share = freqs["C"] / emp_gc if emp_gc > 0 else 0.5
        a_share = freqs["A"] / emp_au if emp_au > 0 else 0.5
        for gc in progress(gc_steps, desc=f"biased gradient {species}", unit="GC"):
            gc_f = float(gc)
            au_f = 1.0 - gc_f
            biased = {
                "C": gc_f * c_share,
                "G": gc_f * (1.0 - c_share),
                "A": au_f * a_share,
                "U": au_f * (1.0 - a_share),
            }
            for _ in range(cfg.gradient_seqs_per_step):
                seq = random_sequence_with_freqs(cfg.sim_seq_length, biased, rng)
                struct, mfe = thermo.fold_mfe(seq)
                length = len(seq)
                gc_pair, au_pair, gu_pair = paired_composition(seq, struct)
                base = _base_composition(seq)
                cond_suffix = " — H-strand" if heavy_strand else ""
                rows.append({
                    "Condition": f"Biased gradient ({species}){cond_suffix}",
                    "Data_Type": "Simulation",
                    "Species": species,
                    "Sequence_GC_Pct": sequence_gc_pct(seq),
                    "GC_Target_Pct": gc_f * 100.0,
                    "Foldedness_Pct": 100.0 * paired_fraction(struct),
                    "Normalized_MFE_per_nt": mfe / length,
                    "Paired_GC_Pct": gc_pair,
                    "Paired_AU_Pct": au_pair,
                    "Paired_GU_Pct": gu_pair,
                    "Pct_A": base["A"],
                    "Pct_C": base["C"],
                    "Pct_G": base["G"],
                    "Pct_U": base["U"],
                })
    return pd.DataFrame(rows)
