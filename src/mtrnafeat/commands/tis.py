"""`mtrnafeat tis` — TIS zoom: −50 / +50 nt around start codon, DMS vs Vienna."""
from __future__ import annotations

from mtrnafeat.analysis import tis
from mtrnafeat.config import Config
from mtrnafeat.io.writers import canonical_csv, tables_csv
from mtrnafeat.viz import tis_plot
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "tis"
    out.mkdir(parents=True, exist_ok=True)
    df = tis.tis_table(cfg, upstream_nt=50, downstream_nt=50)
    canonical_csv(df, out / "tis_dms_vs_mfe.csv")
    tables_csv(df, cfg.outdir, "tis_dms_vs_mfe")
    tis_plot.tis_zoom_panel(df, plot_path(out, "tis_zoom_grid", cfg.plot_format), dpi=cfg.dpi)
    return 0
