"""Per-gene local pair-probability plot, optionally overlaid with the
DMS-derived structure track from the matching .db record.

When DMS overlay data is present, four stacked panels share the x-axis:

* Track 1 — RNAplfold ``P_Paired`` (raw + smoothed). On long transcripts
  the raw layer is faded out so the smoothed line carries the signal.
* Track 2 — DMS-derived paired binary, smoothed at the same rolling
  width as track 1, on the same y-scale.
* Track 3 — Δ = ``mean(P_Paired) − mean(DMS_paired_fraction)`` per
  sliding window, plotted at window centers. Per-window aggregation
  removes the high-frequency oscillation that the per-position
  derivation inherits from the 0/1 DMS step function. Positive = the
  thermodynamic prior predicts more pairing than DMS observes; negative
  = more local accessibility.
* Track 4 — Transcript architecture bar (5'UTR / CDS / 3'UTR), with a
  vertical start-codon line.

A single TIS shading band runs vertically through all four panels so
the eye locks onto the start-codon context immediately. A small
context subtitle (`5'UTR=… · CDS=… · TIS=…`) sits below the title so
the figure can be read on its own.

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
    panel_label,
    style_axis,
)

# A transcript longer than this fades the raw P(paired) layer to barely
# visible — at journal print size the dense per-position spikes form a
# wallpaper of noise that overpowers the smoothed line. The smoothed
# line still tells the story by itself.
_LONG_TRANSCRIPT_NT = 1000
TIS_SHADE_COLOR = "#FFD45A"
TIS_SHADE_ALPHA = 0.18


def _has_dms_overlay(gene_df: pd.DataFrame) -> bool:
    if "DMS_Paired_Smoothed" not in gene_df.columns:
        return False
    return bool(gene_df["DMS_Paired_Smoothed"].notna().any())


def _tis_window_1based(annot: dict, n: int,
                       upstream: int, downstream: int) -> tuple[int, int] | None:
    """Translate the TIS window into 1-based inclusive plot coordinates."""
    cds_start_1 = int(annot["l_utr5"]) + 1
    lo = max(1, cds_start_1 - int(upstream))
    hi = min(n, cds_start_1 + int(downstream))
    if hi <= lo:
        return None
    return lo, hi


def _draw_tis_shade(ax, lo: int, hi: int) -> None:
    ax.axvspan(lo, hi, color=TIS_SHADE_COLOR, alpha=TIS_SHADE_ALPHA, zorder=0)


def _draw_architecture(ax_arch, species: str, gene: str, n: int,
                       annot: dict | None,
                       tis_upstream: int, tis_downstream: int) -> None:
    if annot is None:
        ax_arch.axis("off")
        return
    add_region_track(
        ax_arch,
        l_utr5=int(annot["l_utr5"]),
        l_cds=int(annot["l_cds"]),
        transcript_len=n,
    )
    cds_start_1 = int(annot["l_utr5"]) + 1
    ax_arch.axvline(cds_start_1, color="black", lw=1.0, ls="-", alpha=0.65,
                    zorder=4)


def _context_subtitle(annot: dict | None, window: int, span: int,
                      tis_upstream: int, tis_downstream: int) -> str:
    parts = [f"W={window} nt", f"L={span} nt"]
    if annot is not None:
        l_utr5 = int(annot["l_utr5"])
        l_cds = int(annot["l_cds"])
        l_utr3 = int(annot["l_utr3"])
        parts.append(f"5'UTR={l_utr5} nt")
        parts.append(f"CDS={l_cds} nt")
        parts.append(f"3'UTR={l_utr3} nt")
    parts.append(f"TIS={int(tis_upstream)}/{int(tis_downstream)} nt")
    return "  ·  ".join(parts)


def plot_one_gene(gene_df: pd.DataFrame, out_path: Path,
                  *, smooth_window: int = 25, dpi: int = 300,
                  tis_upstream: int = 50,
                  tis_downstream: int = 50,
                  per_window_df: pd.DataFrame | None = None,
                  scan_window: int | None = None) -> Path:
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

    try:
        annot = annotation_for(species, gene)
    except KeyError:
        annot = None
    tis_window = (
        _tis_window_1based(annot, n, tis_upstream, tis_downstream)
        if annot is not None else None
    )

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
    is_long = n > _LONG_TRANSCRIPT_NT
    raw_alpha = 0.10 if is_long else 0.30
    raw_lw = LINEWIDTH * 0.30 if is_long else LINEWIDTH * 0.45

    # Convention: RNAplfold uses the existing "DMS" palette key (blue) so
    # the legacy plots stay visually consistent. The newly-overlaid
    # *DMS-derived* track is rendered with the RNAstructure brown — the
    # .db dot-bracket was actually produced by RNAstructure upstream, so
    # the color is semantically apt and clearly distinct from the blue.
    rnap_color = PALETTE.get("DMS", "#1F77B4")
    dms_color = PALETTE.get("RNAstructure", "#8C564B")
    delta_color = "#444444"

    if has_dms:
        fig = plt.figure(figsize=(13.5, 8.4))
        gs = fig.add_gridspec(4, 1, height_ratios=[3.0, 3.0, 2.2, 1.0],
                              hspace=0.18)
        ax_p = fig.add_subplot(gs[0])
        ax_dms = fig.add_subplot(gs[1], sharex=ax_p)
        ax_delta = fig.add_subplot(gs[2], sharex=ax_p)
        ax_arch = fig.add_subplot(gs[3], sharex=ax_p)
        panel_label(ax_p, "A")
        panel_label(ax_dms, "B")
        panel_label(ax_delta, "C")

        # Track 1: RNAplfold raw + smoothed
        ax_p.plot(x, raw, color=rnap_color, lw=raw_lw,
                  alpha=raw_alpha, label="P(paired) raw")
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

        # Track 2: DMS paired smoothed
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

        # Track 3: per-window Δ
        # Per-window mean(P_Paired) − mean(DMS_paired_fraction). Aggregated
        # at the same windows used for ``local_probability_per_window.csv``
        # so the plot and the table tell the same story. Falls back to a
        # per-position smoothed Δ when no per_window_df was provided.
        if per_window_df is not None and not per_window_df.empty:
            sub = per_window_df[
                (per_window_df["Species"] == species)
                & (per_window_df["Gene"] == gene)
            ]
            xs = sub["Window_Center_1based"].to_numpy()
            delta = sub["Agreement_Signed_Delta"].to_numpy()
            scan_w = (int(scan_window) if scan_window is not None
                      else int(sub["Window_Size"].median()) if not sub.empty else None)
        else:
            xs = x
            if "Delta_Ppaired_minus_DMS_Smoothed" in gene_df.columns:
                delta = gene_df["Delta_Ppaired_minus_DMS_Smoothed"].to_numpy()
            else:
                delta = smooth - gene_df.get(
                    "DMS_Paired_Smoothed", pd.Series([np.nan] * n)
                ).to_numpy()
            scan_w = None

        finite = np.isfinite(delta)
        mag = float(np.nanmax(np.abs(delta[finite]))) if finite.any() else 0.5
        mag = max(mag, 0.1)
        ax_delta.plot(xs, delta, color=delta_color, lw=LINEWIDTH * 0.95,
                      alpha=0.95)
        ax_delta.fill_between(xs, 0, delta,
                              where=delta >= 0, color=rnap_color,
                              alpha=0.25, interpolate=True,
                              label="↑ RNAplfold > DMS (more pairing)")
        ax_delta.fill_between(xs, 0, delta,
                              where=delta < 0, color=dms_color,
                              alpha=0.25, interpolate=True,
                              label="↓ RNAplfold < DMS (more accessibility)")
        ax_delta.axhline(0, color="black", lw=0.7, ls="-", alpha=0.6)
        ax_delta.set_ylim(-mag * 1.05, mag * 1.05)
        win_tag = f" (per-{scan_w}nt window)" if scan_w else ""
        ax_delta.set_ylabel(f"Δ paired{win_tag}\n(RNAplfold − DMS)",
                            fontsize=LABEL_FONTSIZE)
        ax_delta.tick_params(labelbottom=False)
        ax_delta.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
        ax_delta.set_axisbelow(True)
        legend_outside(ax_delta, position="right", fontsize=9, frameon=False)
        style_axis(ax_delta)

        # Track 4: architecture bar
        _draw_architecture(ax_arch, species, gene, n, annot,
                           int(tis_upstream), int(tis_downstream))
        ax_arch.set_xlim(1, n)
        ax_arch.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
        ax_arch.set_yticks([])
        for sp in ("top", "right", "left"):
            ax_arch.spines[sp].set_visible(False)
        ax_arch.tick_params(axis="y", left=False, labelleft=False)
        style_axis(ax_arch)

        # TIS shading runs through all four panels so the eye locks onto
        # the start-codon context. Drawn at zorder=0, behind the data.
        if tis_window is not None:
            tlo, thi = tis_window
            for ax in (ax_p, ax_dms, ax_delta, ax_arch):
                _draw_tis_shade(ax, tlo, thi)

        ax_p.set_title(
            f"{species} {gene} — RNAplfold local pair probability "
            f"vs DMS-derived structure",
            fontsize=TITLE_FONTSIZE - 1,
            pad=18,
        )
        ax_p.text(
            0.5, 1.02,
            _context_subtitle(annot, window, span, tis_upstream, tis_downstream),
            transform=ax_p.transAxes, ha="center", va="bottom",
            fontsize=LABEL_FONTSIZE - 2, color="#555555",
        )
    else:
        # Legacy 2-panel fallback (no DMS overlay available)
        fig = plt.figure(figsize=(13.5, 5.5))
        gs = fig.add_gridspec(2, 1, height_ratios=[6.0, 1.0], hspace=0.06)
        ax = fig.add_subplot(gs[0])
        ax_arch = fig.add_subplot(gs[1], sharex=ax)

        ax.plot(x, raw, color=rnap_color, lw=raw_lw,
                alpha=raw_alpha, label="P(paired) raw")
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
            f"{species} {gene} — RNAplfold local pair probability",
            fontsize=TITLE_FONTSIZE - 1,
            pad=18,
        )
        ax.text(
            0.5, 1.02,
            _context_subtitle(annot, window, span, tis_upstream, tis_downstream),
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=LABEL_FONTSIZE - 2, color="#555555",
        )
        legend_outside(ax, position="right", fontsize=10, frameon=False)
        style_axis(ax)

        _draw_architecture(ax_arch, species, gene, n, annot,
                           int(tis_upstream), int(tis_downstream))
        ax_arch.set_xlim(1, n)
        ax_arch.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
        ax_arch.set_yticks([])
        for sp in ("top", "right", "left"):
            ax_arch.spines[sp].set_visible(False)
        ax_arch.tick_params(axis="y", left=False, labelleft=False)
        style_axis(ax_arch)

        if tis_window is not None:
            tlo, thi = tis_window
            for axx in (ax, ax_arch):
                _draw_tis_shade(axx, tlo, thi)

    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
