"""`mtrnafeat stats` — per-transcript statistics CSV + boxplot summary."""
from __future__ import annotations

import pandas as pd

from mtrnafeat.analysis.statistics import per_transcript_stats, add_centroid_distances
from mtrnafeat.config import Config
from mtrnafeat.io.writers import canonical_csv, tables_csv
from mtrnafeat.viz import stats_plot
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "stats"
    out.mkdir(parents=True, exist_ok=True)
    frames = []
    for species, fname in cfg.db_files.items():
        df = per_transcript_stats(cfg.data_dir / fname, condition=f"{species} (DMS-guided)")
        df["Species"] = species
        if args and "--centroid" in args:
            df = add_centroid_distances(df, cfg.data_dir / fname)
        frames.append(df)
    out_df = pd.concat(frames, ignore_index=True)
    canonical_csv(out_df, out / "per_transcript_statistics.csv")
    tables_csv(out_df, cfg.outdir, "per_transcript_statistics")
    stats_plot.boxplot_per_species(out_df, plot_path(out, "stats_summary", cfg.plot_format), dpi=cfg.dpi)
    return 0
