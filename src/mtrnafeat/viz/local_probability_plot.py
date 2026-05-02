"""Per-gene local pair-probability plot, optionally overlaid with the
DMS-derived structure track from the matching .db record.

When DMS overlay data is present, four stacked panels share the x-axis:

* Track 1 — RNAplfold ``P_Paired`` raw + smoothed.
* Track 2 — DMS-derived paired binary, smoothed (same rolling-window
  width as track 1) so the two are comparable on the same scale.
* Track 3 — Δ = ``P_Paired_Smoothed − DMS_Paired_Smoothed`` with a
  zero line. Positive = RNAplfold predicts more pairing than DMS;
  negative = RNAplfold predicts more local accessibility.
* Track 4 — Transcript architecture bar (5'UTR / CDS / 3'UTR), with a
  vertical start-codon line and a translucent shaded TIS window.

When DMS overlay data is absent (no matching .db record, or a length
mismatch the analysis layer chose to drop), the plot falls back to the
legacy two-panel layout (P(paired) + architecture bar).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mtrnafeat.constants import PALETTE
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    LINEWIDTH,
    TITLE_FONTSIZE,
    add_region_track,
    apply_theme,
    legend_outside,
    style_axis,
)


def _has_dms_overlay(gene_df: pd.DataFrame) -> bool:
    if "DMS_Paired_Smoothed" not in gene_df.columns:
        return False
    return bool(gene_df["DMS_Paired_Smoothed"].notna().any())


def _draw_architecture(ax_arch, species: str, gene: str, n: int,
                       tis_upstream: int, tis_downstream: int) -> None:
    try:
        annot = annotation_for(species, gene)
    except KeyError:
        ax_arch.axis("off")
        return
    add_region_track(
        ax_arch,
        l_utr5=int(annot["l_utr5"]),
        l_cds=int(annot["l_cds"]),
        transcript_len=n,
    )
    cds_start_1 = int(annot["l_utr5"]) + 1
    ax_arch.axvline(cds_start_1, color="black", lw=1.0, ls="-", alpha=0.55,
                    zorder=4)
    tis_lo = max(1, cds_start_1 - int(tis_upstream))
    tis_hi = min(n, cds_start_1 + int(tis_downstream))
    if tis_hi > tis_lo:
        ax_arch.axvspan(tis_lo, tis_hi, color="#FFD45A", alpha=0.25,
                        zorder=3)


def plot_one_gene(gene_df: pd.DataFrame, out_path: Path,
                  *, smooth_window: int = 25, dpi: int = 300,
                  tis_upstream: int = 50,
                  tis_downstream: int = 50) -> Path:
    apply_theme()
    if gene_df.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "no probability track", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)

    species = str(gene_df["Species"].iloc[0])
    gene = str(gene_df["Gene"].iloc[0])
    window = int(gene_df["Window_nt"].iloc[0])
    span = int(gene_df["Max_BP_Span_nt"].iloc[0])
    n = len(gene_df)
    has_dms = _has_dms_overlay(gene_df)

    x = gene_df["Position_1based"].to_numpy()
    raw = gene_df["P_Paired"].to_numpy()
    if "P_Paired_Smoothed" in gene_df.columns and gene_df["P_Paired_Smoothed"].notna().any():
        smooth = gene_df["P_Paired_Smoothed"].to_numpy()
    else:
        smooth = (
            gene_df["P_Paired"]
            .rolling(int(smooth_window), center=True, min_periods=1)
            .mean()
            .to_numpy()
        )

    # Convention here: RNAplfold uses the existing "DMS" palette key (blue)
    # so the legacy plots stay visually consistent. The newly-overlaid
    # *DMS-derived* track is rendered with the RNAstructure brown — the
    # .db dot-bracket was actually produced by RNAstructure upstream, so
    # the color is semantically apt and clearly distinct from the blue.
    rnap_color = PALETTE.get("DMS", "#1F77B4")
    dms_color = PALETTE.get("RNAstructure", "#8C564B")
    delta_color = "#444444"

    if has_dms:
        fig = plt.figure(figsize=(13.5, 8.0))
        gs = fig.add_gridspec(4, 1, height_ratios=[3.0, 3.0, 2.0, 1.0],
                              hspace=0.18)
        ax_p = fig.add_subplot(gs[0])
        ax_dms = fig.add_subplot(gs[1], sharex=ax_p)
        ax_delta = fig.add_subplot(gs[2], sharex=ax_p)
        ax_arch = fig.add_subplot(gs[3], sharex=ax_p)

        # Track 1: RNAplfold raw + smoothed
        ax_p.plot(x, raw, color=rnap_color, lw=LINEWIDTH * 0.45,
                  alpha=0.35, label="P(paired) raw")
        ax_p.plot(x, smooth, color=rnap_color, lw=LINEWIDTH, alpha=0.95,
                  label=f"P(paired) smoothed (w={smooth_window} nt)")
        ax_p.fill_between(x, 0, smooth, color=rnap_color, alpha=0.15)
        ax_p.set_ylim(-0.02, 1.02)
        ax_p.set_xlim(1, n)
        ax_p.set_ylabel("P(paired)\n(RNAplfold)", fontsize=LABEL_FONTSIZE)
        ax_p.tick_params(labelbottom=False)
        ax_p.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
        ax_p.set_axisbelow(True)
        legend_outside(ax_p, position="right", fontsize=10, frameon=False)
        style_axis(ax_p)

        # Track 2: DMS paired smoothed (and binary as faint stems)
        dms_smoothed = gene_df["DMS_Paired_Smoothed"].to_numpy()
        ax_dms.plot(x, dms_smoothed, color=dms_color, lw=LINEWIDTH,
                    alpha=0.95,
                    label=f"DMS paired fraction (smoothed, w={smooth_window} nt)")
        ax_dms.fill_between(x, 0, dms_smoothed, color=dms_color, alpha=0.15)
        ax_dms.set_ylim(-0.02, 1.02)
        ax_dms.set_ylabel("DMS paired\nfraction", fontsize=LABEL_FONTSIZE)
        ax_dms.tick_params(labelbottom=False)
        ax_dms.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
        ax_dms.set_axisbelow(True)
        legend_outside(ax_dms, position="right", fontsize=10, frameon=False)
        style_axis(ax_dms)

        # Track 3: signed delta with zero line
        delta = gene_df["Delta_Ppaired_minus_DMS_Smoothed"].to_numpy()
        finite = np.isfinite(delta)
        if finite.any():
            mag = float(np.nanmax(np.abs(delta[finite])))
        else:
            mag = 0.5
        mag = max(mag, 0.1)
        ax_delta.plot(x, delta, color=delta_color, lw=LINEWIDTH * 0.85,
                      alpha=0.95)
        ax_delta.fill_between(x, 0, delta,
                              where=delta >= 0, color=rnap_color,
                              alpha=0.20, interpolate=True,
                              label="RNAplfold > DMS (more pairing)")
        ax_delta.fill_between(x, 0, delta,
                              where=delta < 0, color=dms_color,
                              alpha=0.20, interpolate=True,
                              label="RNAplfold < DMS (more openness)")
        ax_delta.axhline(0, color="black", lw=0.7, ls="-", alpha=0.6)
        ax_delta.set_ylim(-mag * 1.05, mag * 1.05)
        ax_delta.set_ylabel("Δ smoothed\n(RNAplfold − DMS)",
                            fontsize=LABEL_FONTSIZE)
        ax_delta.tick_params(labelbottom=False)
        ax_delta.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
        ax_delta.set_axisbelow(True)
        legend_outside(ax_delta, position="right", fontsize=9, frameon=False)
        style_axis(ax_delta)

        # Track 4: architecture bar with TIS shading
        _draw_architecture(ax_arch, species, gene, n,
                           int(tis_upstream), int(tis_downstream))
        ax_arch.set_xlim(1, n)
        ax_arch.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
        ax_arch.set_yticks([])
        for sp in ("top", "right", "left"):
            ax_arch.spines[sp].set_visible(False)
        ax_arch.tick_params(axis="y", left=False, labelleft=False)
        style_axis(ax_arch)

        ax_p.set_title(
            f"{species} {gene} — RNAplfold local pair probability "
            f"vs DMS-derived structure (W={window} nt, L={span} nt)",
            fontsize=TITLE_FONTSIZE - 1,
            pad=12,
        )
    else:
        # Legacy 2-panel fallback (no DMS overlay available)
        fig = plt.figure(figsize=(13.5, 5.5))
        gs = fig.add_gridspec(2, 1, height_ratios=[6.0, 1.0], hspace=0.06)
        ax = fig.add_subplot(gs[0])
        ax_arch = fig.add_subplot(gs[1], sharex=ax)

        ax.plot(x, raw, color=rnap_color, lw=LINEWIDTH * 0.45,
                alpha=0.35, label="P(paired) raw")
        ax.plot(x, smooth, color=rnap_color, lw=LINEWIDTH, alpha=0.95,
                label=f"P(paired) smoothed (w={smooth_window} nt)")
        ax.fill_between(x, 0, smooth, color=rnap_color, alpha=0.15)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlim(1, n)
        ax.set_ylabel("P(paired)\n(RNAplfold)", fontsize=LABEL_FONTSIZE)
        ax.tick_params(labelbottom=False)
        ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
        ax.set_axisbelow(True)
        ax.set_title(
            f"{species} {gene} — RNAplfold local pair probability "
            f"(window={window} nt, max-bp-span={span} nt)",
            fontsize=TITLE_FONTSIZE - 1,
            pad=12,
        )
        legend_outside(ax, position="right", fontsize=10, frameon=False)
        style_axis(ax)

        _draw_architecture(ax_arch, species, gene, n,
                           int(tis_upstream), int(tis_downstream))
        ax_arch.set_xlim(1, n)
        ax_arch.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
        ax_arch.set_yticks([])
        for sp in ("top", "right", "left"):
            ax_arch.spines[sp].set_visible(False)
        ax_arch.tick_params(axis="y", left=False, labelleft=False)
        style_axis(ax_arch)

    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
