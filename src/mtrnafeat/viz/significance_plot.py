"""Refined Panel-B-style visualization for sliding-window dinuc-shuffle z-scores.

Produces:
  - per_window_z_panels: grid of per-(species, gene) panels showing Z_MFE across
    the transcript with peaks (Z < -threshold) marked.
  - peak_overlay: two stacked overlays (Human, Yeast) where every gene's
    Z_MFE trace is drawn against a position normalized to its transcript length,
    so shared structurally significant positions stand out across genes.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from mtrnafeat.constants import PALETTE
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    LINEWIDTH,
    TITLE_FONTSIZE,
    apply_theme,
    shade_regions,
    style_axis,
)


Z_PEAK_THRESHOLD = 2.0
SMOOTH_WINDOW = 5


def _smooth(values: np.ndarray, window: int) -> np.ndarray:
    if len(values) < 2 or window <= 1:
        return values
    s = pd.Series(values).rolling(window, center=True, min_periods=1).mean()
    return s.values


def per_window_z_panels(df: pd.DataFrame, out_path: Path, dpi: int = 300) -> Path:
    """Grid of per-(species, gene) panels: smoothed Z_MFE + peak markers."""
    apply_theme()
    if df is None or df.empty:
        # write an empty placeholder so the file always exists
        fig = plt.figure(figsize=(8, 4))
        plt.text(0.5, 0.5, "no per-window z-scores", ha="center", va="center")
        plt.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)

    keys = sorted(df.groupby(["Species", "Gene"]).groups.keys(),
                   key=lambda sg: (sg[0], sg[1]))
    n = len(keys)
    n_cols = 4 if n >= 8 else max(1, n)
    n_rows = math.ceil(n / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 2.8 * n_rows),
                              squeeze=False)
    axes = axes.flatten()

    for ax, (species, gene) in zip(axes, keys):
        sub = df[(df["Species"] == species) & (df["Gene"] == gene)].copy()
        sub = sub.sort_values("Window_Center_1based")
        x = sub["Window_Center_1based"].values
        z = sub["Z_MFE"].values
        z_smooth = _smooth(z, SMOOTH_WINDOW)

        color = PALETTE.get(species, "#444444")
        ax.plot(x, z, color=color, lw=0.9, alpha=0.35, zorder=2)
        ax.plot(x, z_smooth, color=color, lw=LINEWIDTH * 0.6, alpha=0.95, zorder=3,
                 label="Z_MFE (smoothed)")

        peaks, _ = find_peaks(-z_smooth, height=Z_PEAK_THRESHOLD, distance=8)
        if len(peaks):
            ax.scatter(x[peaks], z_smooth[peaks], color=color, s=55, marker="*",
                        edgecolor="black", linewidths=0.9, zorder=5, label=f"peak (Z < -{Z_PEAK_THRESHOLD})")

        ax.axhline(-Z_PEAK_THRESHOLD, color="black", lw=0.7, ls="--", alpha=0.55, zorder=1)
        ax.axhline(0, color="#888888", lw=0.6, ls=":", alpha=0.6, zorder=1)

        try:
            annot = annotation_for(species, gene)
            shade_regions([ax], l_utr5=annot["l_utr5"], l_cds=annot["l_cds"],
                          transcript_len=annot["l_tr"], label_axis=None)
            ax.set_xlim(1, annot["l_tr"])
        except KeyError:
            pass

        ax.set_title(f"{species} {gene}  (n_peaks={len(peaks)})", fontsize=12, fontweight="bold")
        ax.set_xlabel("Window center (nt)", fontsize=10)
        ax.set_ylabel("Z_MFE", fontsize=10)
        style_axis(ax)

    for ax in axes[len(keys):]:
        ax.set_visible(False)

    fig.suptitle(f"Sliding-window dinuc-shuffle Z_MFE  (peak: Z < -{Z_PEAK_THRESHOLD})",
                  fontsize=TITLE_FONTSIZE, y=1.005)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def peak_overlay(df: pd.DataFrame, out_path: Path, dpi: int = 300) -> Path:
    """Two-row figure (Human top, Yeast bottom): every gene's Z_MFE trace
    overlaid on a normalized position axis so shared peaks across genes stand out.
    """
    apply_theme()
    species_list = ["Human", "Yeast"]
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)

    if df is None or df.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "no per-window z-scores", ha="center", va="center",
                     transform=ax.transAxes)
            ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)

    for ax, species in zip(axes, species_list):
        sub_sp = df[df["Species"] == species]
        if sub_sp.empty:
            ax.text(0.5, 0.5, f"no {species} data", ha="center", va="center",
                     transform=ax.transAxes)
            continue
        genes = sorted(sub_sp["Gene"].unique())
        cmap = plt.get_cmap("tab20" if len(genes) > 10 else "tab10")
        all_peaks_pos: list[float] = []

        for i, gene in enumerate(genes):
            sub = sub_sp[sub_sp["Gene"] == gene].copy().sort_values("Window_Center_1based")
            if len(sub) < 2:
                continue
            x_raw = sub["Window_Center_1based"].values.astype(float)
            x_norm = (x_raw - x_raw.min()) / (x_raw.max() - x_raw.min() + 1e-9)
            z_smooth = _smooth(sub["Z_MFE"].values, SMOOTH_WINDOW)
            color = cmap(i % cmap.N)
            ax.plot(x_norm, z_smooth, color=color, lw=1.2, alpha=0.55, label=gene)
            peaks, _ = find_peaks(-z_smooth, height=Z_PEAK_THRESHOLD, distance=8)
            if len(peaks):
                ax.scatter(x_norm[peaks], z_smooth[peaks], color=color, s=42,
                            edgecolor="black", linewidths=0.7, zorder=5)
                all_peaks_pos.extend(x_norm[peaks].tolist())

        ax.axhline(-Z_PEAK_THRESHOLD, color="black", lw=0.8, ls="--", alpha=0.6)
        ax.axhline(0, color="#888888", lw=0.6, ls=":", alpha=0.6)
        # subtle band at the threshold so the eye latches on to it
        ax.axhspan(-Z_PEAK_THRESHOLD - 0.5, -Z_PEAK_THRESHOLD + 0.5,
                    color="black", alpha=0.04, zorder=0)

        ax.set_title(f"{species} — per-window Z_MFE across {len(genes)} transcripts "
                      f"(normalized position)  ·  {len(all_peaks_pos)} peak(s) total",
                      fontsize=13, fontweight="bold")
        ax.set_ylabel("Z_MFE", fontsize=LABEL_FONTSIZE)
        ax.legend(loc="lower right", fontsize=9, ncol=4, framealpha=0.85)
        style_axis(ax)

    axes[-1].set_xlabel("Normalized transcript position (0 = 5′ end, 1 = 3′ end)",
                        fontsize=LABEL_FONTSIZE)
    fig.suptitle(f"Refined Z_MFE peak overlay (peak: Z < -{Z_PEAK_THRESHOLD})",
                  fontsize=TITLE_FONTSIZE, y=1.005)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
