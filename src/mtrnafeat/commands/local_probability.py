"""`mtrnafeat local-probability` — RNAplfold per-position pair probability,
optionally compared against the DMS-derived dot-bracket from the matching
.db record.

For each gene in ``cfg.target_genes``, runs ViennaRNA's RNAplfold on the
full transcript and emits:

* a per-position pair-probability track (``local_probability_per_position.csv``)
  with optional DMS overlay columns (paired binary, smoothed, signed Δ);
* a sliding-window agreement table (``local_probability_per_window.csv``)
  with mean RNAplfold P(paired) versus DMS paired fraction per window;
* a TIS summary (``local_probability_TIS_summary.csv``) with TIS vs
  CDS-background effect sizes and circular-shift empirical p-values;
* a per-gene plot. When DMS overlay is available the plot has four
  panels (RNAplfold, DMS, Δ, architecture); otherwise it falls back to
  the legacy two-panel layout.

Args (after ``--``):
    --window N              RNAplfold window size (default cfg.rnaplfold_window = 80)
    --max-bp-span N         max base-pair span (default cfg.rnaplfold_max_bp_span = 50)
    --cutoff F              RNAplfold pair-prob cutoff (default cfg.rnaplfold_cutoff = 0.001)
    --smooth N              rolling-window width for smoothed tracks
                            (default cfg.rolling_window)
    --scan-window N         sliding-window size for the per-window agreement
                            table (default cfg.local_probability_scan_window_nt)
    --scan-step N           sliding-window stride (default cfg.local_probability_scan_step_nt)
    --tis-upstream N        TIS upstream context (default cfg.tis_upstream_nt)
    --tis-downstream N      TIS downstream context (default cfg.tis_downstream_nt)
    --tis-n-shuffles N      circular-shift draws for the TIS empirical p
                            (default cfg.tis_n_circular_shifts)
"""
from __future__ import annotations

import pandas as pd

from mtrnafeat.analysis import local_probability
from mtrnafeat.config import Config
from mtrnafeat.constants import file_safe_gene
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.rng import make_rng
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
        elif tok == "--scan-window":
            out["scan_window"] = int(next(it))
        elif tok == "--scan-step":
            out["scan_step"] = int(next(it))
        elif tok == "--tis-upstream":
            out["tis_upstream"] = int(next(it))
        elif tok == "--tis-downstream":
            out["tis_downstream"] = int(next(it))
        elif tok == "--tis-n-shuffles":
            out["tis_n_shuffles"] = int(next(it))
    return out


def run(cfg: Config, args: list[str] | None = None) -> int:
    parsed = _parse(args)
    window = int(parsed.get("window", cfg.rnaplfold_window))
    span = int(parsed.get("max_bp_span", cfg.rnaplfold_max_bp_span))
    cutoff = float(parsed.get("cutoff", cfg.rnaplfold_cutoff))
    smooth = int(parsed.get("smooth", cfg.rolling_window))
    scan_window = int(parsed.get("scan_window", cfg.local_probability_scan_window_nt))
    scan_step = int(parsed.get("scan_step", cfg.local_probability_scan_step_nt))
    tis_upstream = int(parsed.get("tis_upstream", cfg.tis_upstream_nt))
    tis_downstream = int(parsed.get("tis_downstream", cfg.tis_downstream_nt))
    tis_n_shuffles = int(parsed.get("tis_n_shuffles", cfg.tis_n_circular_shifts))

    out = cfg.outdir / "local_probability"
    out.mkdir(parents=True, exist_ok=True)

    df, results = local_probability.scan_all(
        cfg, window=window, max_bp_span=span, cutoff=cutoff,
        smooth=smooth,
        tis_upstream=tis_upstream,
        tis_downstream=tis_downstream,
    )
    if df.empty:
        return 0
    canonical_csv(df, out / "local_probability_per_position.csv")

    # Per-window agreement table — one row per (species, gene, window).
    win_frames = []
    tis_rows = []
    rng = make_rng(cfg.seed + 11)
    for res in results:
        try:
            annot = annotation_for(res.species, res.gene)
        except KeyError:
            annot = None
        win_df = local_probability.per_window_agreement_table(
            res, annot, win=scan_window, step_nt=scan_step,
        )
        if not win_df.empty:
            win_frames.append(win_df)
        if annot is not None:
            tis_rows.append(local_probability.tis_summary_row(
                res, annot,
                upstream=tis_upstream,
                downstream=tis_downstream,
                n_circ_shifts=tis_n_shuffles,
                rng=rng,
            ))

    if win_frames:
        canonical_csv(
            pd.concat(win_frames, ignore_index=True),
            out / "local_probability_per_window.csv",
        )
    if tis_rows:
        canonical_csv(
            pd.DataFrame(tis_rows),
            out / "local_probability_TIS_summary.csv",
        )

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
            tis_upstream=tis_upstream,
            tis_downstream=tis_downstream,
        )
    return 0
