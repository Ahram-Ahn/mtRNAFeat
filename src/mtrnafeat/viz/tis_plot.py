"""TIS zoom visualization: −50 / +50 nt around the start codon.

Two-panel layout (Human left, Yeast right) sharing the y-axis so
cross-species energy comparisons read at a glance.

Visual conventions:

  * x-tick labels are gene names only (single-line, horizontal). The
    available 5'UTR length per gene is shown as a small italic
    annotation directly below each gene's bar group so the metadata
    doesn't bloat the tick label. A small "*" suffix flags the gene
    when the 5'UTR is shorter than the requested upstream window
    (context truncated).
  * Bars whose value is NaN (failed eval / no data) are explicitly
    skipped and the slot is annotated "n/a" above the x-axis baseline.
  * ΔG = 0 (no pairs survived the TIS projection — a *real* signal,
    not missing data) is drawn as a small open diamond at the baseline.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mtrnafeat.viz.samples import ordered_samples
from mtrnafeat.viz.style import LABEL_FONTSIZE, TITLE_FONTSIZE, apply_theme, style_axis

_SPECIES_ORDER = ["Human", "Yeast"]
# Replaced the previous red / blue pairing — the user found that combination
# hard to read. Teal vs amber keeps the two series easily distinguishable
# under all common color-vision deficiencies and reads as paired rather
# than oppositional.
_DMS_COLOR = "#2A9D8F"   # teal
_MFE_COLOR = "#E9A82B"   # amber
_ZERO_MARKER_COLOR = "#222222"


def _bar_with_zero_handling(ax, xs, values, width, color, label, full_ctx):
    """Draw bars; emit a zero-marker diamond when value == 0 and an 'n/a'
    annotation when value is NaN. Truncated-context bars (5'UTR shorter
    than the requested window) are flagged with a "*" suffix on the
    tick label instead of a hatched bar fill, which the user found
    visually noisy."""
    safe_vals = np.where(np.isnan(values), 0.0, values)
    bars = ax.bar(xs, safe_vals, width=width, color=color, label=label,
                   edgecolor="#333333", linewidth=0.7, zorder=3)
    for bar, raw, _full in zip(bars, values, full_ctx, strict=True):
        if np.isnan(raw):
            bar.set_visible(False)
            ax.annotate("n/a", xy=(bar.get_x() + bar.get_width() / 2, 0),
                        xytext=(0, 8), textcoords="offset points",
                        ha="center", va="bottom", fontsize=10, color="#888888",
                        zorder=4)
            continue
        if abs(raw) < 1e-6:
            ax.plot(bar.get_x() + bar.get_width() / 2, 0,
                     marker="D", markersize=6,
                     markerfacecolor="white",
                     markeredgecolor=_ZERO_MARKER_COLOR,
                     markeredgewidth=1.2, zorder=6, clip_on=False)
    return bars


def _draw_one(ax, sub: pd.DataFrame, species: str) -> None:
    if sub.empty:
        ax.text(0.5, 0.5, f"No {species} data.", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="#666666")
        ax.axis("off")
        return
    sub = sub.copy().reset_index(drop=True)
    ind = np.arange(len(sub))
    width = 0.38
    full_ctx = (sub["Has_Full_5UTR_Context"].astype(bool).values
                 if "Has_Full_5UTR_Context" in sub.columns
                 else np.ones(len(sub), dtype=bool))

    _bar_with_zero_handling(
        ax, ind - width / 2, sub["DMS_TIS_Energy"].astype(float).values,
        width, _DMS_COLOR, "DMS-derived", full_ctx,
    )
    _bar_with_zero_handling(
        ax, ind + width / 2, sub["Vienna_TIS_MFE"].astype(float).values,
        width, _MFE_COLOR, "Vienna prediction", full_ctx,
    )

    # Gene name tick labels; * flags truncated 5'UTR context. Labels are
    # drawn diagonally so multi-character mt gene names (ATP8_ATP6,
    # ND4L_ND4) have room and stay legible at journal print size.
    gene_labels = [f"{g}*" if not full else str(g)
                   for g, full in zip(sub["Gene"], full_ctx, strict=True)]
    ax.set_xticks(list(ind))
    ax.set_xticklabels(gene_labels, rotation=45, ha="right",
                        fontsize=LABEL_FONTSIZE + 2, fontweight="bold")
    ax.set_xlim(-0.5, len(sub) - 0.5)
    ax.set_xlabel("Gene  (* = 5'UTR shorter than upstream window)",
                  fontsize=LABEL_FONTSIZE + 1)

    ax.axhline(0, color="black", lw=0.8, zorder=2)
    ax.set_title(species, fontsize=TITLE_FONTSIZE + 1, pad=8, fontweight="bold")
    ax.set_ylabel("ΔG (kcal/mol)", fontsize=LABEL_FONTSIZE + 1)
    ax.tick_params(axis="y", labelsize=LABEL_FONTSIZE)
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
    species_present = ordered_samples(df["Species"].unique())
    n = len(species_present)
    fig, axes = plt.subplots(1, n, figsize=(7.4 * n, 6.4), sharey=True)
    if n == 1:
        axes = [axes]
    for ax, sp in zip(axes, species_present, strict=True):
        _draw_one(ax, df[df["Species"] == sp].sort_values("Gene"), sp)

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=_DMS_COLOR,
                       edgecolor="#333333", linewidth=0.7,
                       label="DMS-derived"),
        plt.Rectangle((0, 0), 1, 1, facecolor=_MFE_COLOR,
                       edgecolor="#333333", linewidth=0.7,
                       label="Vienna prediction"),
        plt.Line2D([0], [0], marker="D", color="none",
                    markerfacecolor="white",
                    markeredgecolor=_ZERO_MARKER_COLOR,
                    markersize=7, markeredgewidth=1.2,
                    label="ΔG = 0 (no pairs)"),
    ]
    fig.legend(handles=handles, loc="upper right", ncol=1,
                frameon=True, framealpha=0.9, fontsize=11)
    fig.suptitle("TIS zoom — DMS vs Vienna ΔG around the start codon (−50 / +50 nt)",
                  fontsize=TITLE_FONTSIZE + 1, y=1.02, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)
