"""Per-gene local pair-probability plot.

Two stacked panels sharing the x-axis:

* Top — RNAplfold ``P_Paired`` track plus a smoothed overlay.
* Bottom — transcript architecture bar (5'UTR / CDS / 3'UTR).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from mtrnafeat.constants import PALETTE
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    LINEWIDTH,
    TITLE_FONTSIZE,
    add_region_track,
    apply_theme,
    style_axis,
)


def plot_one_gene(gene_df: pd.DataFrame, out_path: Path,
                  *, smooth_window: int = 25, dpi: int = 300) -> Path:
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

    fig = plt.figure(figsize=(13.5, 5.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[6.0, 1.0], hspace=0.06)
    ax = fig.add_subplot(gs[0])
    ax_arch = fig.add_subplot(gs[1], sharex=ax)

    x = gene_df["Position_1based"].to_numpy()
    raw = gene_df["P_Paired"].to_numpy()
    smooth = (
        gene_df["P_Paired"]
        .rolling(int(smooth_window), center=True, min_periods=1)
        .mean()
        .to_numpy()
    )

    ax.plot(x, raw, color=PALETTE.get("DMS", "#1F77B4"),
             lw=LINEWIDTH * 0.45, alpha=0.35, label="P(paired) raw")
    ax.plot(x, smooth, color=PALETTE.get("DMS", "#1F77B4"),
             lw=LINEWIDTH, alpha=0.95,
             label=f"P(paired) smoothed (w={smooth_window} nt)")
    ax.fill_between(x, 0, smooth, color=PALETTE.get("DMS", "#1F77B4"),
                     alpha=0.15)

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
    leg = ax.legend(loc="upper right", fontsize=10, frameon=True, framealpha=0.9)
    leg.get_frame().set_edgecolor("#888888")
    leg.get_frame().set_linewidth(0.7)
    style_axis(ax)

    # Architecture bar
    try:
        annot = annotation_for(species, gene)
        add_region_track(
            ax_arch,
            l_utr5=int(annot["l_utr5"]),
            l_cds=int(annot["l_cds"]),
            transcript_len=n,
        )
    except KeyError:
        ax_arch.axis("off")
    ax_arch.set_xlim(1, n)
    ax_arch.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
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
