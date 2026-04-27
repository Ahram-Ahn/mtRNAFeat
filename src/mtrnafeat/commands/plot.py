"""`mtrnafeat plot <kind> --from <dir>` — re-render any plot from cached CSVs."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from mtrnafeat.config import Config
from mtrnafeat.viz import features_plot, landscape_plot
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    if not args:
        print("usage: mtrnafeat plot <kind> [--from <dir>]")
        return 1
    kind = args[0]
    src = cfg.outdir
    if "--from" in args:
        src = Path(args[args.index("--from") + 1])

    if kind == "landscape":
        sim_df = pd.read_csv(src / "landscape" / "specific_conditions.csv")
        grad_df = pd.read_csv(src / "landscape" / "gc_gradient.csv")
        exp_df = pd.read_csv(src / "landscape" / "experimental_overlay.csv")
        landscape_plot.landscape_overlay(
            sim_df, exp_df,
            plot_path(src / "landscape", "landscape_overlay", cfg.plot_format), dpi=cfg.dpi,
        )
        landscape_plot.gradient_curves(
            grad_df,
            plot_path(src / "landscape", "gc_gradient_curves", cfg.plot_format), dpi=cfg.dpi,
        )
        return 0
    if kind == "features":
        motifs = pd.read_csv(src / "features" / "raw_motifs.csv")
        spans = pd.read_csv(src / "features" / "raw_spans.csv")
        features_plot.heatmap_size_ratios(
            motifs, cfg.max_heatmap_size,
            plot_path(src / "features", "heatmap_size_ratios", cfg.plot_format), dpi=cfg.dpi,
        )
        features_plot.phase_space(
            motifs,
            plot_path(src / "features", "phase_space_contour", cfg.plot_format), dpi=cfg.dpi,
        )
        features_plot.span_boxplot(
            spans,
            plot_path(src / "features", "span_boxplot", cfg.plot_format), dpi=cfg.dpi,
        )
        return 0
    print(f"Unknown plot kind: {kind}")
    return 1
