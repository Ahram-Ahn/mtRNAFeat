"""Plots for the structure-deviation stage.

Three figure types::

    plot_one_gene(...)          one PNG/SVG per (species, gene) — four panels
    plot_lollipop(...)          one summary per species — region lollipop
    plot_heatmap(...)           one cross-gene heatmap — bin-mean signed deviation

The colour palette for region classes is fixed across plots so a reader
can recognize a class at a glance regardless of which figure they're
looking at::

    model_high_dms_low   blue       (RNAplfold > DMS;  "DMS-open")
    model_low_dms_high   brown      (DMS > RNAplfold;  "DMS-protected")
    concordant_paired    dark grey  (both high)
    concordant_open      light grey (both low)
    mixed_deviation      purple
    ambiguous            very light grey
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from mtrnafeat.analysis.deviation import DeviationResult
from mtrnafeat.constants import PALETTE
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    LINEWIDTH,
    TITLE_FONTSIZE,
    add_region_track,
    apply_theme,
    legend_outside,
    panel_label,
    style_axis,
)


CLASS_COLORS = {
    "model_high_dms_low": PALETTE.get("DMS", "#1F77B4"),
    "model_low_dms_high": PALETTE.get("RNAstructure", "#8C564B"),
    "concordant_paired": "#444444",
    "concordant_open": "#BFBFBF",
    "mixed_deviation": PALETTE.get("Cotrans", "#9467BD"),
    "ambiguous": "#E0E0E0",
}

CLASS_ORDER = (
    "model_high_dms_low",
    "model_low_dms_high",
    "concordant_paired",
    "concordant_open",
    "mixed_deviation",
    "ambiguous",
)

CLASS_SHORT_LABEL = {
    "model_high_dms_low": "model > DMS (DMS-open)",
    "model_low_dms_high": "DMS > model (DMS-protected)",
    "concordant_paired": "concordant paired",
    "concordant_open": "concordant open",
    "mixed_deviation": "mixed",
    "ambiguous": "ambiguous",
}

TIS_SHADE_COLOR = "#FFD45A"
TIS_SHADE_ALPHA = 0.18


# ──────────────────────── Per-gene figure ────────────────────────


def _draw_class_legend(ax, classes_present: list[str]) -> None:
    handles = []
    labels = []
    for cls in CLASS_ORDER:
        if cls in classes_present:
            handles.append(Rectangle((0, 0), 1, 1,
                                     color=CLASS_COLORS[cls], alpha=0.55))
            labels.append(CLASS_SHORT_LABEL[cls])
    if not handles:
        return
    ax.legend(handles, labels,
              loc="upper center", bbox_to_anchor=(0.5, -0.55),
              ncol=min(3, len(handles)),
              frameon=False, fontsize=9, borderaxespad=0.0)


def _gene_subtitle(annot: dict | None, result: DeviationResult,
                   tis_up: int, tis_down: int) -> str:
    parts = [
        f"W={result.rnaplfold_window} nt",
        f"L={result.rnaplfold_max_bp_span} nt",
    ]
    if annot is not None:
        parts.append(f"5'UTR={int(annot['l_utr5'])} nt")
        parts.append(f"CDS={int(annot['l_cds'])} nt")
        parts.append(f"3'UTR={int(annot['l_utr3'])} nt")
    parts.append(f"TIS=−{int(tis_up)}/+{int(tis_down)} nt")
    return "  ·  ".join(parts)


def plot_one_gene(result: DeviationResult,
                  regions_df: pd.DataFrame,
                  out_path: Path,
                  *, cfg, dpi: int = 300,
                  per_window_df: pd.DataFrame | None = None,
                  scan_window: int | None = None) -> Path:
    apply_theme()
    n = len(result.sequence)
    species, gene = result.species, result.gene
    try:
        annot = annotation_for(species, gene)
    except KeyError:
        annot = None

    rnap_color = PALETTE.get("DMS", "#1F77B4")
    dms_color = PALETTE.get("RNAstructure", "#8C564B")

    fig = plt.figure(figsize=(13.5, 8.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[3.0, 3.0, 2.8, 1.0],
                          hspace=0.18)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1], sharex=ax_a)
    ax_c = fig.add_subplot(gs[2], sharex=ax_a)
    ax_d = fig.add_subplot(gs[3], sharex=ax_a)
    panel_label(ax_a, "A")
    panel_label(ax_b, "B")
    panel_label(ax_c, "C")
    panel_label(ax_d, "D")

    x = np.arange(1, n + 1)

    # Panel A — RNAplfold P_paired
    ax_a.plot(x, result.p_model_smooth, color=rnap_color, lw=LINEWIDTH,
              alpha=0.95, label="P(paired) model (smoothed)")
    ax_a.fill_between(x, 0, result.p_model_smooth, color=rnap_color, alpha=0.15)
    ax_a.set_ylim(-0.02, 1.02)
    ax_a.set_xlim(1, n)
    ax_a.set_ylabel("Model\nP(paired)", fontsize=LABEL_FONTSIZE)
    ax_a.tick_params(labelbottom=False)
    ax_a.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    ax_a.set_axisbelow(True)
    style_axis(ax_a)

    # Panel B — DMS paired fraction
    ax_b.plot(x, result.p_dms_smooth, color=dms_color, lw=LINEWIDTH,
              alpha=0.95, label="DMS paired fraction (smoothed)")
    ax_b.fill_between(x, 0, result.p_dms_smooth, color=dms_color, alpha=0.15)
    ax_b.set_ylim(-0.02, 1.02)
    ax_b.set_ylabel("DMS\npaired", fontsize=LABEL_FONTSIZE)
    ax_b.tick_params(labelbottom=False)
    ax_b.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    ax_b.set_axisbelow(True)
    style_axis(ax_b)

    # Panel C — signed deviation. We render at the per-window scale that
    # ``local-probability`` uses (window=cfg.local_probability_scan_window_nt,
    # step=…) so the two stages tell visually identical Δ stories. Region
    # calling itself still runs on the per-position 25-nt rolling track
    # in the analysis layer; the rectangles below carry that information.
    use_per_window = (per_window_df is not None and not per_window_df.empty)
    if use_per_window:
        sub = per_window_df[
            (per_window_df["Species"] == species)
            & (per_window_df["Gene"] == gene)
        ]
        xs_dev = sub["Window_Center_1based"].to_numpy()
        dev = sub["Agreement_Signed_Delta"].to_numpy()
        scan_w = (int(scan_window) if scan_window is not None
                  else int(sub["Window_Size"].median()) if not sub.empty else None)
    else:
        xs_dev = x
        dev = result.deviation_smooth
        scan_w = None

    finite = np.isfinite(dev)
    mag = float(np.nanmax(np.abs(dev[finite]))) if finite.any() else 0.5
    mag = max(mag, 0.3)
    ax_c.plot(xs_dev, dev, color="#222222", lw=LINEWIDTH * 0.95)
    ax_c.fill_between(xs_dev, 0, dev, where=dev >= 0, color=rnap_color,
                      alpha=0.18, interpolate=True)
    ax_c.fill_between(xs_dev, 0, dev, where=dev < 0, color=dms_color,
                      alpha=0.18, interpolate=True)
    ax_c.axhline(0, color="black", lw=0.7, alpha=0.6)
    # Threshold lines are only meaningful at the per-position scale that
    # region calling uses; on the per-window line they'd visually undercut
    # the called rectangles, so we omit them when the per-window frame is
    # supplied.
    if not use_per_window:
        thr = float(cfg.structure_deviation_threshold)
        ax_c.axhline(thr, color="black", lw=0.6, ls=":", alpha=0.6)
        ax_c.axhline(-thr, color="black", lw=0.6, ls=":", alpha=0.6)
    # Region rectangles colored by class. Visual encoding is already
    # carried by (a) these tinted bands, (b) the colored region blocks
    # on the architecture bar (panel D), and (c) the bottom legend.
    # Inline class-name annotations stack on top of each other on
    # gene-dense transcripts and add no information — omitted on
    # purpose.
    classes_present: set[str] = set()
    for _, r in regions_df.iterrows():
        cls = r["Region_Class"]
        classes_present.add(cls)
        color = CLASS_COLORS.get(cls, "#888888")
        rect = Rectangle(
            (int(r["Start_1based"]), -mag * 1.05),
            int(r["End_1based"]) - int(r["Start_1based"]),
            mag * 2.10,
            color=color, alpha=0.22, zorder=0,
        )
        ax_c.add_patch(rect)
    ax_c.set_ylim(-mag * 1.05, mag * 1.05)
    win_tag = f" (per-{scan_w}nt window)" if scan_w else ""
    ax_c.set_ylabel(f"Signed deviation{win_tag}\n(model − DMS)",
                    fontsize=LABEL_FONTSIZE)
    ax_c.tick_params(labelbottom=False)
    ax_c.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    ax_c.set_axisbelow(True)
    # Compact direction key on the right of panel C
    ax_c.text(1.005, 0.95, "↑ model > DMS",
              transform=ax_c.transAxes, fontsize=9, color=rnap_color,
              ha="left", va="top")
    ax_c.text(1.005, 0.05, "↓ DMS > model",
              transform=ax_c.transAxes, fontsize=9, color=dms_color,
              ha="left", va="bottom")
    style_axis(ax_c)

    # Panel D — gene architecture + region blocks colored by class
    if annot is not None:
        add_region_track(
            ax_d,
            l_utr5=int(annot["l_utr5"]),
            l_cds=int(annot["l_cds"]),
            transcript_len=n,
        )
        cds_start_1 = int(annot["l_utr5"]) + 1
        ax_d.axvline(cds_start_1, color="black", lw=1.0, ls="-",
                     alpha=0.65, zorder=4)
    # Region blocks drawn just above the architecture bar
    for _, r in regions_df.iterrows():
        cls = r["Region_Class"]
        color = CLASS_COLORS.get(cls, "#888888")
        rect = Rectangle(
            (int(r["Start_1based"]), 0.55),
            int(r["End_1based"]) - int(r["Start_1based"]),
            0.30,
            color=color, alpha=0.85, zorder=5,
        )
        ax_d.add_patch(rect)
    ax_d.set_xlim(1, n)
    ax_d.set_ylim(-0.05, 1.05)
    ax_d.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
    ax_d.set_yticks([])
    for sp in ("top", "right", "left"):
        ax_d.spines[sp].set_visible(False)
    ax_d.tick_params(axis="y", left=False, labelleft=False)
    style_axis(ax_d)

    # TIS shading vertically through every panel
    if annot is not None:
        cds_start_1 = int(annot["l_utr5"]) + 1
        tis_up = int(cfg.structure_deviation_tis_upstream)
        tis_down = int(cfg.structure_deviation_tis_downstream)
        tis_lo = max(1, cds_start_1 - tis_up)
        tis_hi = min(n, cds_start_1 + tis_down)
        if tis_hi > tis_lo:
            for ax in (ax_a, ax_b, ax_c, ax_d):
                ax.axvspan(tis_lo, tis_hi, color=TIS_SHADE_COLOR,
                           alpha=TIS_SHADE_ALPHA, zorder=0)

    # Title + subtitle
    ax_a.set_title(
        f"{species} {gene} — DMS / model structure-deviation regions",
        fontsize=TITLE_FONTSIZE - 1, pad=18,
    )
    ax_a.text(
        0.5, 1.02,
        _gene_subtitle(annot, result,
                       cfg.structure_deviation_tis_upstream,
                       cfg.structure_deviation_tis_downstream),
        transform=ax_a.transAxes, ha="center", va="bottom",
        fontsize=LABEL_FONTSIZE - 2, color="#555555",
    )

    # Single legend (bottom of the figure) with the classes that are
    # actually called for this gene.
    _draw_class_legend(ax_d, sorted(classes_present, key=lambda c: CLASS_ORDER.index(c)
                                    if c in CLASS_ORDER else 999))

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


# ──────────────────────── Lollipop ────────────────────────


def plot_lollipop(regions_df: pd.DataFrame, out_path: Path,
                  *, species: str, cfg, dpi: int = 300) -> Path | None:
    """One row per gene; lollipops at each region's midpoint, height =
    mean deviation, dot size ∝ length, color = class."""
    apply_theme()
    sub = regions_df[regions_df["Species"] == species].copy()
    if sub.empty:
        return None
    genes = sorted(sub["Gene"].unique().tolist())
    n_genes = len(genes)
    fig_h = max(2.5, 0.55 * n_genes + 1.4)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    gene_y = {g: i for i, g in enumerate(genes)}

    classes_present: set[str] = set()
    for _, r in sub.iterrows():
        y = gene_y[r["Gene"]]
        cls = r["Region_Class"]
        classes_present.add(cls)
        color = CLASS_COLORS.get(cls, "#888888")
        height = float(r["Mean_Deviation"])
        size = max(20.0, min(320.0, float(r["Length_nt"]) * 1.0))
        # Stem from y to y + height (height in deviation units, but we
        # want the stem to read as "how strong" without overlapping rows).
        # Trick: scale stem length to a fixed fraction of row height.
        stem_top = y + np.sign(height) * 0.42 * (abs(height) / 0.75)
        ax.plot([float(r["Midpoint_1based"]), float(r["Midpoint_1based"])],
                [y, stem_top], color=color, lw=1.4, alpha=0.85, zorder=3)
        ax.scatter([float(r["Midpoint_1based"])], [stem_top],
                   s=size, color=color, edgecolor="white", linewidths=0.6,
                   alpha=0.9, zorder=4)

    ax.set_yticks(list(gene_y.values()))
    ax.set_yticklabels(genes, fontsize=LABEL_FONTSIZE - 1)
    ax.set_ylim(-0.5, n_genes - 0.5)
    ax.invert_yaxis()
    ax.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
    ax.set_title(
        f"{species} — structure-deviation regions per gene "
        f"(stem ∝ |mean Δ|, dot size ∝ length)",
        fontsize=TITLE_FONTSIZE - 1, pad=10,
    )
    ax.grid(True, axis="x", linestyle="--", linewidth=0.5, alpha=0.35)
    ax.set_axisbelow(True)
    style_axis(ax)

    # Legend outside on the right
    handles = []
    labels = []
    for cls in CLASS_ORDER:
        if cls in classes_present:
            handles.append(plt.Line2D([0], [0], marker="o", color="white",
                                      markerfacecolor=CLASS_COLORS[cls],
                                      markeredgecolor="white",
                                      markersize=10, lw=0))
            labels.append(CLASS_SHORT_LABEL[cls])
    if handles:
        ax.legend(handles, labels, loc="center left",
                  bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=9)

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


# ──────────────────────── Heatmap ────────────────────────


_BIN_ORDER = ("5_end", "TIS", "early_CDS", "mid_CDS",
              "late_CDS", "stop_proximal", "3_end")


def plot_heatmap(matrix_df: pd.DataFrame, out_path: Path,
                 *, value: str = "Mean_Deviation",
                 cfg=None, dpi: int = 300) -> Path | None:
    """Cross-gene heatmap. Rows = (Species, Gene); columns = region bin.

    ``value`` is the column to read from ``matrix_df``; defaults to
    signed deviation (the publication-headline metric).
    """
    apply_theme()
    if matrix_df.empty:
        return None
    df = matrix_df.copy()
    df["Row_Label"] = df["Species"] + " " + df["Gene"]
    pivot = df.pivot_table(
        index="Row_Label", columns="Region_Bin", values=value, aggfunc="mean"
    )
    cols = [b for b in _BIN_ORDER if b in pivot.columns]
    pivot = pivot[cols]
    # Sort rows: yeast first then human, alphabetically within species
    species_rank = {"Yeast": 0, "Human": 1}
    pivot = pivot.reindex(sorted(
        pivot.index,
        key=lambda r: (species_rank.get(r.split(" ", 1)[0], 99),
                       r.split(" ", 1)[1] if " " in r else r),
    ))

    n_rows, n_cols = pivot.shape
    fig_h = max(3.5, 0.45 * n_rows + 1.6)
    fig_w = max(7.0, 1.1 * n_cols + 4.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Symmetric colormap centered on zero for signed deviation; otherwise
    # use a sequential map.
    if value == "Mean_Deviation":
        vmax = float(np.nanmax(np.abs(pivot.values))) if pivot.size else 0.5
        vmax = max(vmax, 0.1)
        im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax)
        cbar_label = "Mean Δ paired (model − DMS)"
    else:
        im = ax.imshow(pivot.values, aspect="auto", cmap="viridis")
        cbar_label = value.replace("_", " ")
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(cols, rotation=30, ha="right",
                       fontsize=LABEL_FONTSIZE - 1)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(pivot.index, fontsize=LABEL_FONTSIZE - 1)
    ax.set_title(
        f"Cross-gene structure deviation per architectural bin "
        f"({cbar_label})",
        fontsize=TITLE_FONTSIZE - 1, pad=10,
    )
    # Annotate each cell with the value
    for i in range(n_rows):
        for j in range(n_cols):
            v = pivot.values[i, j]
            if not np.isfinite(v):
                continue
            # Choose text color for legibility on RdBu_r
            txt_color = "white" if (value == "Mean_Deviation" and abs(v) > 0.55 * vmax) else "black"
            ax.text(j, i, f"{v:+.2f}" if value == "Mean_Deviation" else f"{v:.2f}",
                    ha="center", va="center", fontsize=8, color=txt_color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(cbar_label, fontsize=LABEL_FONTSIZE - 1)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)
