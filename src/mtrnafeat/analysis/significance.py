"""Dinucleotide-shuffle z-score per gene and per sliding window.

Workman & Krogh 1999 framework: a sequence's MFE is "significantly stable" if
it is more negative than the MFE of dinucleotide-shuffled controls. Same
dinucleotide composition isolates structural signal from compositional bias.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mtrnafeat.analysis.window import sliding_intervals
from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core import thermo
from mtrnafeat.core.shuffle import dinuc_shuffle
from mtrnafeat.io.db_parser import parse_db
from mtrnafeat.progress import progress, step
from mtrnafeat.rng import make_rng


def gene_zscore(seq: str, n_shuffles: int, rng: np.random.Generator,
                 label: str | None = None) -> dict:
    _, observed = thermo.fold_mfe(seq)
    null = np.empty(n_shuffles, dtype=float)
    bar = progress(range(n_shuffles), desc=label, unit="shuffle", leave=False) if label else range(n_shuffles)
    for k in bar:
        sh = dinuc_shuffle(seq, rng)
        _, mfe = thermo.fold_mfe(sh)
        null[k] = mfe
    mu = float(np.mean(null))
    sigma = float(np.std(null, ddof=1))
    z = (observed - mu) / sigma if sigma > 0 else float("nan")
    p = float(np.mean(null <= observed))  # one-sided P(null ≤ observed)
    return {
        "MFE_observed": observed,
        "Null_Mean": mu,
        "Null_Std": sigma,
        "Z_MFE": z,
        "P_Empirical": p,
        "N_Shuffles": int(n_shuffles),
    }


def per_gene_significance(cfg: Config) -> pd.DataFrame:
    rng = make_rng(cfg.seed + 7)
    step(f"running significance per gene (n_shuffles={cfg.n_shuffles})")
    rows = []
    all_recs = [(species, rec)
                for species, fname in cfg.db_files.items()
                for rec in parse_db(cfg.data_dir / fname)]
    for species, rec in progress(all_recs, desc="significance (genes)", unit="gene"):
        res = gene_zscore(rec.sequence, cfg.n_shuffles, rng,
                           label=f"{species} {rec.gene}")
        rows.append({"Species": species, "Gene": canonical_gene(rec.gene),
                     "Length": len(rec.sequence), **res})
    return pd.DataFrame(rows)


def per_window_significance(cfg: Config) -> pd.DataFrame:
    """Sliding-window z-score: ScanFold-style per-window null."""
    rng = make_rng(cfg.seed + 8)
    step(f"running significance per window (n_shuffles={cfg.n_shuffles})")
    rows = []
    for species, fname in cfg.db_files.items():
        path = cfg.data_dir / fname
        for rec in progress(list(parse_db(path)),
                              desc=f"{species} per-window", unit="gene"):
            n = len(rec.sequence)
            intervals = sliding_intervals(n, cfg.window_nt, cfg.step_nt)
            for s, e in progress(intervals, desc=f"{species} {rec.gene} windows",
                                  unit="win", leave=False):
                seq_w = rec.sequence[s:e]
                if len(seq_w) < 30:
                    continue
                res = gene_zscore(seq_w, cfg.n_shuffles, rng)
                rows.append({
                    "Species": species,
                    "Gene": canonical_gene(rec.gene),
                    "Window_Start_1based": s + 1,
                    "Window_End_1based": e,
                    "Window_Center_1based": s + 1 + (len(seq_w) - 1) / 2.0,
                    **res,
                })
    return pd.DataFrame(rows)
