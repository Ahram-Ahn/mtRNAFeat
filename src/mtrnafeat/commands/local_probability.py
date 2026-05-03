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
from mtrnafeat.core.stats import bh_fdr
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.rng import make_rng
from mtrnafeat.viz import local_probability_plot
from mtrnafeat.viz.style import plot_path

_TIS_PVAL_COLUMNS: tuple[tuple[str, str], ...] = (
    ("Empirical_P_TIS_Low_Ppaired_vs_CDS_Windows",
     "Q_Value_TIS_Low_Ppaired_vs_CDS_Windows"),
    ("Empirical_P_TIS_Low_DMS_Pairfrac_vs_CDS_Windows",
     "Q_Value_TIS_Low_DMS_Pairfrac_vs_CDS_Windows"),
)


def _add_bh_q_columns(df: pd.DataFrame, group_col: str | None = None) -> pd.DataFrame:
    """Add Q_Value_* columns next to each Empirical_P_* column.

    BH-FDR is applied within ``group_col`` if provided (e.g. per
    sensitivity-sweep window so each TIS context width is its own
    hypothesis family), otherwise across the whole table (the primary
    TIS summary, where one row per gene defines the family).
    """
    if df.empty:
        for _, q_col in _TIS_PVAL_COLUMNS:
            df[q_col] = pd.Series(dtype=float)
        return df
    for p_col, q_col in _TIS_PVAL_COLUMNS:
        if p_col not in df.columns:
            continue
        if group_col is None:
            df[q_col] = bh_fdr(df[p_col].to_numpy(dtype=float))
        else:
            df[q_col] = (
                df.groupby(group_col, sort=False)[p_col]
                .transform(lambda s: bh_fdr(s.to_numpy(dtype=float)))
            )
    return df


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
        else:
            raise SystemExit(
                f"local-probability: unknown flag {tok!r}. "
                "See module docstring for the supported flag list."
            )
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
    sens_rows: list[dict] = []
    rng = make_rng(cfg.seed + 11)
    sweep_pairs = tuple((int(u), int(d)) for u, d in cfg.tis_window_sweep_pairs)
    primary_pair = (int(tis_upstream), int(tis_downstream))
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
            primary_row = local_probability.tis_summary_row(
                res, annot,
                upstream=tis_upstream,
                downstream=tis_downstream,
                n_circ_shifts=tis_n_shuffles,
                rng=rng,
            )
            tis_rows.append(primary_row)
            # Sensitivity sweep: same machinery at multiple TIS context
            # widths so reviewers can judge robustness rather than rely on
            # one fixed window. Reuses the primary row when the (u, d)
            # pair matches so we don't double-pay the circular-shift cost.
            for u, d in sweep_pairs:
                if (u, d) == primary_pair:
                    sens_rows.append({
                        "TIS_Window_Tag": f"{u}_{d}",
                        **primary_row,
                    })
                    continue
                row = local_probability.tis_summary_row(
                    res, annot,
                    upstream=u, downstream=d,
                    n_circ_shifts=tis_n_shuffles,
                    rng=rng,
                )
                sens_rows.append({"TIS_Window_Tag": f"{u}_{d}", **row})

    per_window_all = (pd.concat(win_frames, ignore_index=True)
                      if win_frames else pd.DataFrame())
    if not per_window_all.empty:
        canonical_csv(
            per_window_all,
            out / "local_probability_per_window.csv",
        )
    if tis_rows:
        tis_df = _add_bh_q_columns(pd.DataFrame(tis_rows), group_col=None)
        canonical_csv(tis_df, out / "local_probability_TIS_summary.csv")
    if sens_rows:
        # Each TIS context width (e.g. ±30, ±50, ±100, ±200, ±500) is its
        # own hypothesis family — apply BH within window tag, not across.
        sens_df = _add_bh_q_columns(
            pd.DataFrame(sens_rows), group_col="TIS_Window_Tag",
        )
        canonical_csv(sens_df, out / "local_probability_TIS_sensitivity.csv")

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
            per_window_df=per_window_all,
            scan_window=scan_window,
        )
    return 0
