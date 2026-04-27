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
from mtrnafeat.core.shuffle import random_gc_sequence, random_sequence_with_freqs
from mtrnafeat.core.structure import paired_fraction
from mtrnafeat.io.db_parser import parse_db
from mtrnafeat.progress import progress, step
from mtrnafeat.rng import make_rng


def compute_empirical_freqs(db_path) -> dict[str, float]:
    """Aggregate (A, U, G, C) frequency over every record in a .db file."""
    counts: Counter[str] = Counter()
    for rec in parse_db(db_path):
        counts.update(rec.sequence.upper())
    total = sum(counts[b] for b in "AUGC")
    if total == 0:
        return {"A": 0.25, "U": 0.25, "G": 0.25, "C": 0.25}
    return {b: counts[b] / total for b in "AUGC"}


def species_freqs_for_pipeline(cfg: Config) -> dict[str, dict[str, float]]:
    """Return {species: {A,U,G,C}} for every species declared in cfg.db_files,
    using either cfg.sim_freqs_per_species (if provided) or empirical .db data."""
    out: dict[str, dict[str, float]] = {}
    given = getattr(cfg, "sim_freqs_per_species", None) or {}
    for species, fname in cfg.db_files.items():
        if species in given:
            out[species] = {b: float(given[species].get(b, 0.0)) for b in "AUGC"}
        else:
            out[species] = compute_empirical_freqs(cfg.data_dir / fname)
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
            gc_pair, au_pair, gu_pair = paired_composition(rec.sequence, rec.structure)
            rows.append({
                "Condition": f"Exp: {species} DMS",
                "Data_Type": "Experimental",
                "Gene": rec.gene, "Species": species, "Length": length,
                "MFE": rec.mfe,
                "Normalized_MFE_per_nt": rec.mfe / length if length else 0.0,
                "Foldedness_Pct": 100.0 * paired_fraction(rec.structure),
                "Sequence_GC_Pct": sequence_gc_pct(rec.sequence),
                "Paired_GC_Pct": gc_pair,
                "Paired_AU_Pct": au_pair,
                "Paired_GU_Pct": gu_pair,
            })
    return pd.DataFrame(rows)
