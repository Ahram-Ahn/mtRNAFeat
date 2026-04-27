"""`mtrnafeat significance` — dinucleotide-shuffle z-score per gene + per window."""
from __future__ import annotations

from mtrnafeat.analysis import significance
from mtrnafeat.config import Config
from mtrnafeat.io.writers import canonical_csv


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "significance"
    out.mkdir(parents=True, exist_ok=True)
    canonical_csv(significance.per_gene_significance(cfg), out / "z_per_gene.csv")
    if args and "--per-window" in args:
        canonical_csv(significance.per_window_significance(cfg), out / "z_per_window.csv")
    return 0
