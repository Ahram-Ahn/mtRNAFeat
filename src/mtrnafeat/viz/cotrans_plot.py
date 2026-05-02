"""Per-gene cotranscriptional / sliding-window signal plots.

One PNG per (species, gene): four stacked panels sharing the x-axis.

* MFE / nt — energetic stability of the prefix or window
* Ensemble diversity — Vienna's mean BP distance, an uncertainty proxy
* Paired fraction — fraction of nt paired in the MFE structure
* Concordance — z-scored deltas of all three, with peak markers and
  background shading where two or more signals peak together.

Background tinting on the concordance panel: green where MFE and diversity
*both* drop sharply (a "structure-commits" event), red where both rise (a
"structure-loosens" event). Mirrors the legacy ``signal_analysis.py``
shading idiom.
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
    apply_theme,
    style_axis,
)


def _peaks_in(z: np.ndarray, threshold: float, sign: int, distance: int = 5) -> np.ndarray:
    from scipy.signal import find_peaks
    if sign < 0:
        peaks, _ = find_peaks(-z, height=threshold, distance=distance)
    else:
        peaks, _ = find_peaks(z, height=threshold, distance=distance)
    return peaks


def _shade_concordance(ax, x: np.ndarray, z_mfe: np.ndarray, z_div: np.ndarray) -> None:
    """Shade x-ranges where MFE and diversity both drop (or both rise) sharply."""
    if x.size < 2:
        return
    drop_band = (z_mfe < -0.6) & (z_div < -0.6)
    rise_band = (z_mfe > 0.6) & (z_div > 0.6)
    y0, y1 = ax.get_ylim()
    ax.fill_between(x, y0, y1, where=drop_band, color="forestgreen",
                     alpha=0.15, zorder=0, step="mid")
    ax.fill_between(x, y0, y1, where=rise_band, color="crimson",
                     alpha=0.10, zorder=0, step="mid")


def plot_one_gene(gene_df: pd.DataFrame, out_path: Path,
                  *, z_threshold: float = 2.0, dpi: int = 300) -> Path:
    """Four-panel per-gene plot for the cotranscriptional / sliding signals."""
    apply_theme()
    if gene_df.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "no signal", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)

    species = str(gene_df["Species"].iloc[0])
    gene = str(gene_df["Gene"].iloc[0])
    mode = str(gene_df["Mode"].iloc[0])

    # X axis: end-of-window position (matches "nascent transcript length"
    # interpretation for prefix mode and is intuitive for sliding mode too).
    x = gene_df["Window_End_1based"].to_numpy(dtype=float)

    fig, axes = plt.subplots(4, 1, figsize=(13.5, 10.0), sharex=True,
                              gridspec_kw={"height_ratios": [1, 1, 1, 1.25]})
    species_color = PALETTE.get(species, "#444444")

    # Panel 1 — MFE / nt
    ax = axes[0]
    ax.plot(x, gene_df["MFE_per_nt"], color=species_color, lw=LINEWIDTH * 0.85)
    ax.set_ylabel("MFE / nt\n(kcal/mol)", fontsize=LABEL_FONTSIZE - 1)
    ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    ax.set_axisbelow(True)
    style_axis(ax)

    # Panel 2 — ensemble diversity
    ax = axes[1]
    ax.plot(x, gene_df["Diversity"], color=PALETTE.get("Diversity", "#000000"),
            lw=LINEWIDTH * 0.85)
    ax.set_ylabel("Ensemble\ndiversity (nt)", fontsize=LABEL_FONTSIZE - 1)
    ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    ax.set_axisbelow(True)
    style_axis(ax)

    # Panel 3 — paired fraction
    ax = axes[2]
    ax.plot(x, gene_df["Paired_Fraction"], color=PALETTE.get("DMS", "#1F77B4"),
            lw=LINEWIDTH * 0.85)
    ax.set_ylabel("Paired\nfraction", fontsize=LABEL_FONTSIZE - 1)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    ax.set_axisbelow(True)
    style_axis(ax)

    # Panel 4 — z-score concordance
    ax = axes[3]
    have_z = ("Z_Delta_MFE_per_nt_Smooth" in gene_df.columns
              and "Z_Delta_Diversity_Smooth" in gene_df.columns
              and "Z_Delta_Paired_Fraction_Smooth" in gene_df.columns)
    if have_z:
        z_mfe = gene_df["Z_Delta_MFE_per_nt_Smooth"].to_numpy()
        z_div = gene_df["Z_Delta_Diversity_Smooth"].to_numpy()
        z_pf = gene_df["Z_Delta_Paired_Fraction_Smooth"].to_numpy()
        ax.plot(x, z_mfe, color=species_color, lw=1.4, alpha=0.9,
                 label="Z(ΔMFE/nt)")
        ax.plot(x, z_div, color=PALETTE.get("Diversity", "#000000"),
                 lw=1.4, alpha=0.9, label="Z(ΔDiversity)")
        ax.plot(x, z_pf, color=PALETTE.get("DMS", "#1F77B4"),
                 lw=1.4, alpha=0.9, label="Z(ΔPaired_Fraction)")
        ax.axhline(0, color="#888888", lw=0.6, ls=":", alpha=0.7)
        ax.axhline(-z_threshold, color="black", lw=0.7, ls="--", alpha=0.55)
        ax.axhline(z_threshold, color="black", lw=0.7, ls="--", alpha=0.55)

        # Peak markers
        for arr, sign, color, label in (
            (z_mfe, -1, species_color, "ΔMFE peak"),
            (z_div, -1, PALETTE.get("Diversity", "#000000"), "ΔDiv peak"),
            (z_pf, 1, PALETTE.get("DMS", "#1F77B4"), "ΔPF peak"),
        ):
            peaks = _peaks_in(arr, threshold=z_threshold, sign=sign)
            if peaks.size:
                ax.scatter(x[peaks], arr[peaks], color=color, s=70, marker="*",
                            edgecolor="black", linewidths=0.9, zorder=6,
                            label=f"{label} (n={peaks.size})")

        # Concordance shading (computed AFTER ylim is finalized for the y0/y1 read)
        ax.relim(); ax.autoscale_view()
        _shade_concordance(ax, x, z_mfe, z_div)

        # Custom placement: the cotrans bottom panel already carries the
        # shared X-axis label of the 4-panel stack, so the standard
        # `legend_outside(..., position="bottom")` collides with it. Push
        # the legend below the X-axis label.
        ax.legend(
            loc="upper center", bbox_to_anchor=(0.5, -0.36),
            fontsize=9, ncol=4, frameon=False, borderaxespad=0.0,
        )
        ax.set_ylabel("Z-score\n(smoothed Δ)", fontsize=LABEL_FONTSIZE - 1)
    style_axis(ax)

    # X-axis label
    if mode == "prefix":
        axes[-1].set_xlabel("Nascent transcript length (nt)", fontsize=LABEL_FONTSIZE)
    else:
        axes[-1].set_xlabel("Window end position (nt)", fontsize=LABEL_FONTSIZE)

    # Title with sample counts
    n_pts = len(gene_df)
    n_peaks = 0
    if have_z:
        n_peaks = sum(
            int((-arr if sign < 0 else arr).max() >= z_threshold)
            for arr, sign in (
                (gene_df["Z_Delta_MFE_per_nt_Smooth"].to_numpy(), -1),
                (gene_df["Z_Delta_Diversity_Smooth"].to_numpy(), -1),
                (gene_df["Z_Delta_Paired_Fraction_Smooth"].to_numpy(), 1),
            )
            if arr.size
        )
    fig.suptitle(
        f"{species} {gene} — cotranscriptional signals "
        f"(mode={mode}, n={n_pts} samples, peak |Z|≥{z_threshold:.1f})",
        fontsize=TITLE_FONTSIZE - 1,
        fontweight="bold",
        y=0.995,
    )

    # Try to draw region boundaries from the annotation table
    try:
        annot = annotation_for(species, gene)
        cds_start = annot["l_utr5"] + 0.5
        cds_end = annot["l_utr5"] + annot["l_cds"] + 0.5
        for ax in axes:
            ax.axvline(cds_start, color="#888888", ls="--", lw=0.7, alpha=0.55,
                       zorder=0)
            ax.axvline(cds_end, color="#888888", ls="--", lw=0.7, alpha=0.55,
                       zorder=0)
    except KeyError:
        pass

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
