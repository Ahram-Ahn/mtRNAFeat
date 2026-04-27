"""`mtrnafeat significance` — dinucleotide-shuffle z-score per gene + per window.

Args (after `--`):
    --per-window     additionally compute per-window z-scores and emit the
                     refined Panel-B-style plots. Slow (≈n_shuffles × n_windows
                     folds per gene); `mtrnafeat run-all` passes this by default.
"""
from __future__ import annotations

from mtrnafeat.analysis import significance
from mtrnafeat.config import Config
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import significance_plot
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "significance"
    out.mkdir(parents=True, exist_ok=True)
    canonical_csv(significance.per_gene_significance(cfg), out / "z_per_gene.csv")

    if args and "--per-window" in args:
        df_win = significance.per_window_significance(cfg)
        canonical_csv(df_win, out / "z_per_window.csv")
        significance_plot.per_window_z_panels(
            df_win, plot_path(out, "z_per_window_panels", cfg.plot_format), dpi=cfg.dpi)
        significance_plot.peak_overlay(
            df_win, plot_path(out, "z_peak_overlay", cfg.plot_format), dpi=cfg.dpi)
    return 0
