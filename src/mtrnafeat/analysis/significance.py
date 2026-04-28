"""Dinucleotide-shuffle z-score per gene.

Workman & Krogh 1999 framework: a sequence's MFE is "significantly stable" if
it is more negative than the MFE of dinucleotide-shuffled controls. Same
dinucleotide composition isolates structural signal from compositional bias.

The per-window equivalent now lives in ``analysis.cotrans`` (richer signals,
per-gene plots).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

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


def _records_for_targets(cfg: Config) -> list[tuple[str, "object"]]:
    """Yield (species, DbRecord) for every record whose canonical gene name
    is in ``cfg.target_genes``. If ``target_genes`` is empty/None, fall back
    to every record in every .db file."""
    targets = {canonical_gene(g) for g in (cfg.target_genes or ())}
    out: list[tuple[str, object]] = []
    for species, fname in cfg.db_files.items():
        for rec in parse_db(cfg.data_dir / fname):
            if not targets or canonical_gene(rec.gene) in targets:
                out.append((species, rec))
    return out


def per_gene_significance(cfg: Config) -> pd.DataFrame:
    rng = make_rng(cfg.seed + 7)
    step(f"running significance per gene (n_shuffles={cfg.n_shuffles})")
    rows = []
    all_recs = _records_for_targets(cfg)
    for species, rec in progress(all_recs, desc="significance (genes)", unit="gene"):
        res = gene_zscore(rec.sequence, cfg.n_shuffles, rng,
                           label=f"{species} {rec.gene}")
        rows.append({"Species": species, "Gene": canonical_gene(rec.gene),
                     "Length": len(rec.sequence), **res})
    return pd.DataFrame(rows)


# Per-window dinuc-shuffle z-scoring used to live here. It's been retired
# in favour of `analysis.cotrans` — same goal (find structurally interesting
# positions across each transcript) but with richer Vienna-derived signals
# (ensemble diversity, paired fraction, MFE/nt) and per-gene plots.
