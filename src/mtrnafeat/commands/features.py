"""`mtrnafeat features` — element decomposition + heatmap + phase space + box."""
from __future__ import annotations

import pandas as pd

from mtrnafeat.analysis import features
from mtrnafeat.config import Config
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import features_plot
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "features"
    out.mkdir(parents=True, exist_ok=True)

    dms_motifs, dms_spans = features.features_dms(cfg)
    sim_motifs, sim_spans = features.features_simulated(cfg)
    motifs = pd.concat([dms_motifs, sim_motifs], ignore_index=True)
    spans = pd.concat([dms_spans, sim_spans], ignore_index=True)
    region = features.region_stratified_elements(cfg)

    canonical_csv(motifs, out / "raw_motifs.csv")
    canonical_csv(spans, out / "raw_spans.csv")
    canonical_csv(region, out / "region_stratified_elements.csv")

    features_plot.heatmap_size_ratios(motifs, cfg.max_heatmap_size,
                                       plot_path(out, "heatmap_size_ratios", cfg.plot_format), dpi=cfg.dpi)
    features_plot.phase_space(motifs, plot_path(out, "phase_space_contour", cfg.plot_format), dpi=cfg.dpi)
    features_plot.span_boxplot(spans, plot_path(out, "span_boxplot", cfg.plot_format), dpi=cfg.dpi)
    return 0
