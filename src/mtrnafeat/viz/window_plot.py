"""Whole-transcript paired-fraction plot for the `window` command.

Two stacked panels share the x-axis (transcript position):

* **Top** — rolling paired-fraction track (centered window, default 25 nt).
  Three lines: DMS-guided structure (black), Vienna full (blue), and the
  engine-span fold (orange when Vienna, brown when RNAstructure). ΔG values
  go in the legend so the figure reads on its own.
* **Bottom** — transcript architecture bar (5'UTR / CDS / 3'UTR), with
  region labels placed inside the wider bands and via leader lines for
  any band <5% of transcript length.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from mtrnafeat.constants import PALETTE
from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    LEGEND_FONTSIZE,
    LINEWIDTH,
    TITLE_FONTSIZE,
    add_region_track,
    apply_theme,
    style_axis,
)


def _engine_label(engine: str) -> str:
    return "Vienna" if engine == "vienna" else "RNAstructure"


def _engine_color(engine: str) -> str:
    return PALETTE.get("Vienna_span", "#FF7F0E") if engine == "vienna" \
        else PALETTE.get("RNAstructure", "#8C564B")


def plot_transcript_pairing(res, pos_df: pd.DataFrame, out_path: Path,
                            rolling_window: int, dpi: int = 300) -> Path:
    """Two-panel figure: rolling paired-fraction track + transcript architecture bar."""
    apply_theme()

    eng = _engine_label(res.engine)
    span = res.max_bp_span
    engine_color = _engine_color(res.engine)

    dms_col = "DMS_RollingPairedFrac"
    vfull_col = "ViennaFull_RollingPairedFrac"
    espan_col = f"{eng}Span{span}_RollingPairedFrac"

    n = len(res.sequence)

    fig = plt.figure(figsize=(15.5, 7.0))
    gs = fig.add_gridspec(
        2, 1, height_ratios=[6.5, 1.0], hspace=0.05,
    )
    ax = fig.add_subplot(gs[0])
    ax_arch = fig.add_subplot(gs[1], sharex=ax)

    x = pos_df["Position_1based"].to_numpy()

    ax.plot(
        x, pos_df[dms_col],
        color=PALETTE.get("DMS", "#1F77B4"),
        lw=LINEWIDTH,
        label=f"DMS-guided  (ΔG = {res.dms_recalc_mfe:,.1f} kcal/mol)",
    )
    ax.plot(
        x, pos_df[vfull_col],
        color=PALETTE.get("Vienna", "#D62728"),
        lw=LINEWIDTH * 0.9,
        alpha=0.95,
        label=f"Vienna full  (ΔG = {res.vienna_full_mfe:,.1f} kcal/mol)",
    )
    ax.plot(
        x, pos_df[espan_col],
        color=engine_color,
        lw=LINEWIDTH * 0.9,
        alpha=0.95,
        label=f"{eng} max-bp-span = {span} nt  "
              f"(ΔG = {res.engine_span_mfe_native:,.1f} kcal/mol)",
    )

    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(1, n)
    ax.set_ylabel(
        f"Paired fraction\n({rolling_window}-nt rolling window, 1-nt step)",
        fontsize=LABEL_FONTSIZE,
    )
    title = (
        f"{res.species} {res.gene}: local pairing profile across the full transcript"
    )
    ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=14)
    ax.tick_params(labelbottom=False)
    ax.grid(True, axis="y", linestyle="--", linewidth=0.6, alpha=0.35)
    ax.set_axisbelow(True)
    leg = ax.legend(
        loc="upper right", frameon=True, framealpha=0.92,
        fontsize=LEGEND_FONTSIZE, borderpad=0.6, handlelength=2.2,
    )
    leg.get_frame().set_edgecolor("#888888")
    leg.get_frame().set_linewidth(0.7)
    style_axis(ax)

    add_region_track(
        ax_arch,
        l_utr5=int(res.annot["l_utr5"]),
        l_cds=int(res.annot["l_cds"]),
        transcript_len=n,
    )
    ax_arch.set_xlim(1, n)
    ax_arch.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
    # Re-show x-axis on the architecture panel even though add_region_track
    # turned the axis off.
    ax_arch.set_axis_on()
    ax_arch.set_yticks([])
    ax_arch.spines["top"].set_visible(False)
    ax_arch.spines["right"].set_visible(False)
    ax_arch.spines["left"].set_visible(False)
    ax_arch.tick_params(axis="y", left=False, labelleft=False)
    style_axis(ax_arch)

    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
