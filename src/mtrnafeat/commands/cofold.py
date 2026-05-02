"""`mtrnafeat cofold` — CoFold parameter sweep against DMS ΔG.

For each (species, gene), folds the wild-type CDS at every (alpha, tau)
combo on the configured grid (cfg.cofold_alpha_sweep × cfg.cofold_tau_sweep)
and reports how close each combo gets to the DMS-evaluated ΔG. Also runs
a per-window CoFold-vs-DMS correlation analysis.

Args (after `--`):
    --no-window-corr     skip the per-window correlation pass (faster).
"""
from __future__ import annotations

from mtrnafeat.analysis import cofold_sweep
from mtrnafeat.config import Config
from mtrnafeat.io.writers import canonical_csv, tables_csv
from mtrnafeat.viz import cofold_plot


def _parse(args: list[str] | None) -> dict:
    out: dict = {"do_window_corr": True}
    if not args:
        return out
    for tok in args:
        if tok == "--no-window-corr":
            out["do_window_corr"] = False
    return out


def run(cfg: Config, args: list[str] | None = None) -> int:
    parsed = _parse(args)
    out = cfg.outdir / "cofold"
    out.mkdir(parents=True, exist_ok=True)

    full, win = cofold_sweep.run_cofold_sweep(cfg, do_window_corr=parsed["do_window_corr"])
    if full.empty:
        print("[mtrnafeat] cofold: no genes matched. Skipping.")
        return 0
    canonical_csv(full, out / "cofold_grid.csv")
    best = cofold_sweep.best_per_gene(full)
    if not best.empty:
        canonical_csv(best, out / "cofold_best_per_gene.csv")
        tables_csv(best, cfg.outdir, "cofold_best_per_gene")
    if not win.empty:
        canonical_csv(win, out / "cofold_per_window_corr.csv")

    cofold_plot.gap_strip_panels(full, out, cfg.plot_format, dpi=cfg.dpi)
    if not win.empty:
        cofold_plot.per_window_corr_curves(win, out, cfg.plot_format, dpi=cfg.dpi)
    return 0
