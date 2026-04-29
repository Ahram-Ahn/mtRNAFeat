"""TIS zoom visualization: −50 / +50 nt around the start codon.

Publication-quality two-panel layout: Human (left) and Yeast (right),
sharing the y-axis so cross-species energy comparisons read at a glance.

Visual conventions designed to keep ΔG = 0 from being confused with
"missing data":

  * A **distinct ΔG = 0 marker** (a small open diamond) is drawn at the
    bar position when the energy is exactly zero — this happens when the
    DMS dot-bracket has no pairs surviving the projection onto the TIS
    window, which is a *real* signal ("the in-vivo structure here is
    unfolded") and must be visually different from a missing bar.
  * Bars whose value is NaN (failed eval / no data) are explicitly
    skipped and the slot is annotated "n/a" above the x-axis baseline.
  * Bars sitting next to a ΔG = 0 marker are drawn at full opacity so
    the contrast between "zero, real" and "missing" is unambiguous.

Hatched bars still flag genes whose 5'UTR is shorter than `upstream_nt`
(upstream context truncated). x-tick labels show the actual usable
upstream length per gene as `GENE (5'UTR=N)`.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mtrnafeat.viz.style import LABEL_FONTSIZE, TITLE_FONTSIZE, apply_theme, style_axis

_SPECIES_ORDER = ["Human", "Yeast"]
_DMS_COLOR = "#1F4E79"
_MFE_COLOR = "#C0392B"
_ZERO_MARKER_COLOR = "#222222"


def _bar_with_zero_handling(ax, xs, values, width, color, label, full_ctx):
    """Draw bars; emit a zero-marker diamond when value == 0 and an 'n/a'
    annotation when value is NaN. Returns the BarContainer (without the
    NaN slots, which would otherwise paint as zero-height rectangles)."""
    safe_vals = np.where(np.isnan(values), 0.0, values)
    bars = ax.bar(xs, safe_vals, width=width, color=color, label=label,
                   edgecolor="#333333", linewidth=0.7, zorder=3)
    for bar, raw, full in zip(bars, values, full_ctx):
        if np.isnan(raw):
            bar.set_visible(False)
            ax.annotate("n/a", xy=(bar.get_x() + bar.get_width() / 2, 0),
                        xytext=(0, 8), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8, color="#888888",
                        zorder=4)
            continue
        if not full:
            bar.set_hatch("///")
            bar.set_edgecolor("white")
        if abs(raw) < 1e-6:
            # Open diamond at the zero baseline so "ΔG = 0" reads
            # differently from a missing bar.
            ax.plot(bar.get_x() + bar.get_width() / 2, 0,
                     marker="D", markersize=7,
                     markerfacecolor="white",
                     markeredgecolor=_ZERO_MARKER_COLOR,
                     markeredgewidth=1.4, zorder=6, clip_on=False)
            ax.annotate("0", xy=(bar.get_x() + bar.get_width() / 2, 0),
                         xytext=(0, 12), textcoords="offset points",
                         ha="center", va="bottom", fontsize=8,
                         color=_ZERO_MARKER_COLOR, zorder=6)
    return bars


def _draw_one(ax, sub: pd.DataFrame, species: str) -> None:
    if sub.empty:
        ax.text(0.5, 0.5, f"No {species} data.", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#666666")
        ax.axis("off")
        return
    sub = sub.copy().reset_index(drop=True)
    ind = np.arange(len(sub))
    width = 0.4
    full_ctx = (sub["Has_Full_5UTR_Context"].astype(bool).values
                 if "Has_Full_5UTR_Context" in sub.columns
                 else np.ones(len(sub), dtype=bool))

    _bar_with_zero_handling(
        ax, ind - width / 2, sub["DMS_TIS_Energy"].astype(float).values,
        width, _DMS_COLOR, "DMS-projected ΔG (Vienna eval)", full_ctx,
    )
    _bar_with_zero_handling(
        ax, ind + width / 2, sub["Vienna_TIS_MFE"].astype(float).values,
        width, _MFE_COLOR, "Vienna MFE", full_ctx,
    )

    if "L_5UTR_in_window" in sub.columns:
        labels = [f"{g}\n(5'UTR={int(u)} nt)" for g, u in zip(sub["Gene"], sub["L_5UTR_in_window"])]
    else:
        labels = list(sub["Gene"])
    ax.set_xticks(list(ind))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=10)
    ax.axhline(0, color="black", lw=0.8, zorder=2)
    ax.set_title(species, fontsize=TITLE_FONTSIZE - 1, pad=8, fontweight="bold")
    ax.set_ylabel("ΔG (kcal/mol)", fontsize=LABEL_FONTSIZE)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.6, alpha=0.45)
    ax.set_axisbelow(True)
    style_axis(ax)


def tis_zoom_panel(df_tis: pd.DataFrame, out_path: Path, dpi: int = 300) -> Path:
    apply_theme()
    if df_tis.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No TIS data.", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)
    df = df_tis.copy()
    # Keep rows even if one of the two energies is NaN — the bar-drawer
    # surfaces missing values explicitly with an "n/a" annotation rather
    # than silently dropping the gene.
    df = df.dropna(subset=["DMS_TIS_Energy", "Vienna_TIS_MFE"], how="all")
    species_present = [s for s in _SPECIES_ORDER if s in df["Species"].unique()]
    if not species_present:
        species_present = sorted(df["Species"].unique())
    n = len(species_present)
    fig, axes = plt.subplots(1, n, figsize=(7.4 * n, 5.8), sharey=True)
    if n == 1:
        axes = [axes]
    for ax, sp in zip(axes, species_present):
        _draw_one(ax, df[df["Species"] == sp].sort_values("Gene"), sp)

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=_DMS_COLOR,
                       edgecolor="#333333", linewidth=0.7,
                       label="DMS-projected ΔG (Vienna eval on .db structure)"),
        plt.Rectangle((0, 0), 1, 1, facecolor=_MFE_COLOR,
                       edgecolor="#333333", linewidth=0.7, label="Vienna MFE"),
        plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="#333333",
                       linewidth=0.7, hatch="///",
                       label="5'UTR < 50 nt (upstream context truncated)"),
        plt.Line2D([0], [0], marker="D", color="none",
                    markerfacecolor="white",
                    markeredgecolor=_ZERO_MARKER_COLOR,
                    markersize=8, markeredgewidth=1.4,
                    label="ΔG = 0 (no projected pairs — real, not missing)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2,
                bbox_to_anchor=(0.5, -0.06), frameon=True, framealpha=0.95,
                fontsize=10)
    fig.suptitle("TIS zoom — DMS-projected ΔG vs Vienna MFE around the start codon (−50 / +50 nt)",
                  fontsize=TITLE_FONTSIZE - 1, y=1.02, fontweight="bold")
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)
