"""Whole-transcript paired-fraction plot for the `window` command.

Publication-quality two-panel layout. Panels share the x-axis (transcript
position) and the legend lives **outside** the data axis.

* **Top** — rolling paired-fraction track (centered window, default
  25 nt). Two lines: the DMS-guided structure and a thermodynamic
  prediction at the configured max-bp-span (default 300 nt). The
  full-length Vienna fold without a max-bp-span cap was previously
  shown as a third trace; it added clutter without information for
  long mRNAs (max_bp_span = 300 already captures every realistic
  contact at this scale) and has been removed. ΔG values appear next
  to each label so the figure reads on its own.
* **Bottom** — transcript architecture bar (5'UTR / CDS / 3'UTR) with
  consistent region labels (inside wide bands, leader-lined when
  narrow).
* **Right margin** — legend, drawn in its own figure-coordinate slot so
  no plot real-estate is overlapped at any window size.
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
    """Two-panel figure: rolling paired-fraction trace + transcript architecture bar.

    The legend is placed outside the trace axis (top-right, in figure
    coordinates) so the three pairing-fraction traces are never occluded
    by it on long-transcript or peak-rich genes.
    """
    apply_theme()

    eng = _engine_label(res.engine)
    span = res.max_bp_span
    engine_color = _engine_color(res.engine)

    dms_col = "DMS_RollingPairedFrac"
    espan_col = f"{eng}Span{span}_RollingPairedFrac"

    n = len(res.sequence)

    # Wider canvas + an extra figure-side margin for the right-hand legend.
    fig = plt.figure(figsize=(17.0, 6.4))
    gs = fig.add_gridspec(
        2, 1, height_ratios=[6.5, 1.0], hspace=0.06,
        left=0.06, right=0.78, top=0.88, bottom=0.13,
    )
    ax = fig.add_subplot(gs[0])
    ax_arch = fig.add_subplot(gs[1], sharex=ax)

    x = pos_df["Position_1based"].to_numpy()

    ax.plot(
        x, pos_df[dms_col],
        color=PALETTE.get("DMS", "#1F77B4"),
        lw=LINEWIDTH,
        label=f"DMS-derived  (ΔG = {res.dms_recalc_mfe:,.1f} kcal/mol)",
    )
    ax.plot(
        x, pos_df[espan_col],
        color=engine_color,
        lw=LINEWIDTH * 0.9,
        alpha=0.95,
        label=f"{eng} prediction  (ΔG = {res.engine_span_mfe_native:,.1f} kcal/mol)",
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
    ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=12, fontweight="bold")
    ax.tick_params(labelbottom=False)
    ax.grid(True, axis="y", linestyle="--", linewidth=0.6, alpha=0.35)
    ax.set_axisbelow(True)
    # Put the legend OUTSIDE the data axis (figure top-right). bbox is in
    # figure coords because we constrained gs.right to 0.78.
    leg = fig.legend(
        loc="upper left", bbox_to_anchor=(0.79, 0.88),
        frameon=False,
        fontsize=LEGEND_FONTSIZE, borderpad=0.7, handlelength=2.4,
        labelspacing=0.9, title=f"DMS vs {eng}",
        title_fontsize=LEGEND_FONTSIZE,
    )
    style_axis(ax)

    add_region_track(
        ax_arch,
        l_utr5=int(res.annot["l_utr5"]),
        l_cds=int(res.annot["l_cds"]),
        transcript_len=n,
    )
    ax_arch.set_xlim(1, n)
    ax_arch.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
    ax_arch.set_axis_on()
    ax_arch.set_yticks([])
    ax_arch.spines["top"].set_visible(False)
    ax_arch.spines["right"].set_visible(False)
    ax_arch.spines["left"].set_visible(False)
    ax_arch.tick_params(axis="y", left=False, labelleft=False)
    style_axis(ax_arch)

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)
