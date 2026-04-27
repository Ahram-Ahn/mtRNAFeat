"""Sliding-window full-transcript plot with region shading.

One figure per (species, gene). The x-axis is always the full transcript
(1..transcript_len). UTR vs CDS vs 3'UTR are conveyed as subtle
background tints + dashed boundary lines + region labels above the top
panel — short 5'UTRs (common in human mt-mRNA) get a leader line
pointing into the narrow band so the label is still legible.

The CDS-only zoom plot was retired — it was redundant with the
full-transcript figure for every gene we care about, and it forced the
user to compare two figures side by side.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from mtrnafeat.constants import PALETTE
from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    LINEWIDTH,
    TITLE_FONTSIZE,
    apply_theme,
    shade_regions,
    style_axis,
)


def plot_full_transcript(window_df, l_utr5: int, l_cds: int, transcript_len: int,
                          species: str, gene: str, max_bp_span: int | None,
                          out_path: Path, dpi: int = 300) -> Path:
    """Two-panel figure: ΔG and paired-fraction tracks across the full transcript."""
    apply_theme()
    fig = plt.figure(figsize=(15.5, 8.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.18)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    x = window_df["Window_Center_1based"].values

    span_label = f" (max_bp_span={max_bp_span})" if max_bp_span else ""

    ax1.plot(x, window_df["DMS_Window_Energy"], color=PALETTE.get("DMS", "#1F77B4"),
             lw=LINEWIDTH, label="DMS-projected ΔG")
    ax1.plot(x, window_df["Vienna_Window_MFE"], color=PALETTE.get("Vienna", "#D62728"),
             lw=LINEWIDTH, label="Vienna local MFE")
    ax1.set_ylabel("Window ΔG (kcal/mol)", fontsize=LABEL_FONTSIZE)
    ax1.set_title(f"{species} {gene} — sliding-window scan{span_label}",
                   fontsize=TITLE_FONTSIZE, pad=18)
    ax1.legend(loc="best", fontsize=11)
    style_axis(ax1)

    ax2.plot(x, window_df["DMS_Window_PairedFraction"], color=PALETTE.get("DMS", "#1F77B4"),
             lw=LINEWIDTH, label="DMS paired fraction")
    ax2.plot(x, window_df["Vienna_Window_PairedFraction"], color=PALETTE.get("Vienna", "#D62728"),
             lw=LINEWIDTH, label="Vienna paired fraction")
    ax2.set_ylabel("Window paired fraction", fontsize=LABEL_FONTSIZE)
    ax2.set_xlabel("Transcript position (window center, nt)", fontsize=LABEL_FONTSIZE)
    ax2.set_ylim(0, 1)
    ax2.legend(loc="best", fontsize=11)
    style_axis(ax2)

    shade_regions([ax1, ax2], l_utr5=l_utr5, l_cds=l_cds, transcript_len=transcript_len,
                    label_axis=ax1)
    ax1.set_xlim(1, transcript_len)
    ax2.set_xlim(1, transcript_len)

    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
