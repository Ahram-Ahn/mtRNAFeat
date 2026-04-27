"""`mtrnafeat substitution` — synonymous-recoding thermodynamic permutation test.

Promoted from `legacy/base_substitution/03.*`. Three null pools per gene
(flat-GC, positional-GC, synonymous) compared against the wild-type ΔG
under CoFold. Emits both the long-format distribution CSV and a per-gene
summary table; figures are KDE small-multiples + a Z-heatmap.

Args (after `--`):
    --n INT          override cfg.substitution_n_simulations (default: 200)
    --max-nt INT     override cfg.substitution_max_nt (default: 300)
"""
from __future__ import annotations

from dataclasses import replace

from mtrnafeat.analysis import substitution_thermo
from mtrnafeat.config import Config
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import substitution_plot


def _parse(args: list[str] | None) -> dict:
    out: dict = {}
    if not args:
        return out
    it = iter(args)
    for tok in it:
        if tok == "--n":
            out["substitution_n_simulations"] = int(next(it))
        elif tok == "--max-nt":
            out["substitution_max_nt"] = int(next(it))
    return out


def run(cfg: Config, args: list[str] | None = None) -> int:
    overrides = _parse(args)
    if overrides:
        cfg = replace(cfg, **overrides)
    out = cfg.outdir / "substitution"
    out.mkdir(parents=True, exist_ok=True)
    tables = cfg.outdir / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    dist, summary = substitution_thermo.run_substitution_thermo(cfg)
    if dist.empty:
        print("[mtrnafeat] substitution: no genes had codon-complete sequences ≥60 nt; skipping.")
        return 0
    canonical_csv(dist, out / "substitution_thermo_distribution.csv")
    canonical_csv(summary, tables / "substitution_thermo_summary.csv")
    substitution_plot.kde_panels(dist, out, cfg.plot_format, dpi=cfg.dpi)
    substitution_plot.z_heatmap(summary, out, cfg.plot_format, dpi=cfg.dpi)
    return 0
