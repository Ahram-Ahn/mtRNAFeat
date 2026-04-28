"""`mtrnafeat local-probability` — RNAplfold per-position pair probability.

For each gene in ``cfg.target_genes``, runs ViennaRNA's RNAplfold on the
full transcript and emits a per-position pair-probability track.

Args (after ``--``):
    --window N         RNAplfold window size (default cfg.rnaplfold_window = 80)
    --max-bp-span N    max base-pair span (default cfg.rnaplfold_max_bp_span = 50)
    --cutoff F         RNAplfold pair-prob cutoff (default cfg.rnaplfold_cutoff = 0.001)
    --smooth N         rolling-window width for the smoothed plot overlay
                       (default cfg.rolling_window)
"""
from __future__ import annotations

from mtrnafeat.analysis import local_probability
from mtrnafeat.config import Config
from mtrnafeat.constants import file_safe_gene
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import local_probability_plot
from mtrnafeat.viz.style import plot_path


def _parse(args: list[str] | None) -> dict:
    out: dict = {}
    if not args:
        return out
    it = iter(args)
    for tok in it:
        if tok == "--window":
            out["window"] = int(next(it))
        elif tok == "--max-bp-span":
            out["max_bp_span"] = int(next(it))
        elif tok == "--cutoff":
            out["cutoff"] = float(next(it))
        elif tok == "--smooth":
            out["smooth"] = int(next(it))
    return out


def run(cfg: Config, args: list[str] | None = None) -> int:
    parsed = _parse(args)
    window = int(parsed.get("window", cfg.rnaplfold_window))
    span = int(parsed.get("max_bp_span", cfg.rnaplfold_max_bp_span))
    cutoff = float(parsed.get("cutoff", cfg.rnaplfold_cutoff))
    smooth = int(parsed.get("smooth", cfg.rolling_window))

    out = cfg.outdir / "local_probability"
    out.mkdir(parents=True, exist_ok=True)

    df, results = local_probability.scan_all(
        cfg, window=window, max_bp_span=span, cutoff=cutoff
    )
    if df.empty:
        return 0
    canonical_csv(df, out / "local_probability_per_position.csv")

    for res in results:
        gene_df = df[(df["Species"] == res.species) & (df["Gene"] == res.gene)]
        local_probability_plot.plot_one_gene(
            gene_df.reset_index(drop=True),
            out_path=plot_path(
                out,
                f"local_probability_{res.species}_{file_safe_gene(res.gene)}",
                cfg.plot_format,
            ),
            smooth_window=smooth,
            dpi=cfg.dpi,
        )
    return 0
