"""Yeast↔human COX1 comparative plots: substitution bar chart + directional flux."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from mtrnafeat.viz.style import LABEL_FONTSIZE, TICK_FONTSIZE, TITLE_FONTSIZE, apply_theme


def plot_directional_flux(flux_df, out_path: Path, dpi: int = 300) -> Path:
    """12 (from→to) × 3 (codon position) heatmap of bias score.

    Rows sorted by |mean bias| across positions (most biased on top).
    BH-significant cells outlined in black (q<0.05).
    """
    apply_theme()
    if flux_df.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No directional flux data.", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)
    df = flux_df.copy()
    df["Direction"] = df["From"] + "→" + df["To"]
    pivot = df.pivot_table(index="Direction", columns="Position", values="Bias_Score", aggfunc="mean")
    qpivot = df.pivot_table(index="Direction", columns="Position", values="BH_q", aggfunc="mean")

    # Sort rows by absolute mean bias across all positions (most biased on top)
    row_order = pivot.abs().mean(axis=1).sort_values(ascending=False).index
    pivot = pivot.loc[row_order]
    qpivot = qpivot.reindex(row_order)

    fig, ax = plt.subplots(figsize=(7.5, 8))
    sns.heatmap(pivot, ax=ax, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                annot_kws={"size": 10}, linewidths=0.4,
                cbar_kws={"label": "Bias score (obs − exp)"})

    # Outline BH-significant cells
    for i, direction in enumerate(pivot.index):
        for j, pos in enumerate(pivot.columns):
            try:
                q = qpivot.loc[direction, pos]
            except KeyError:
                q = 1.0
            if q < 0.05:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                                            edgecolor="black", lw=2.0))

    ax.set_title(
        "Yeast→Human substitution flux (BH-significant outlined, q<0.05)\n"
        "Positive = enriched;  negative = depleted vs null expectation",
        fontsize=TITLE_FONTSIZE - 2, pad=10,
    )
    ax.set_xlabel("Codon position  (1 = first,  2 = second,  3 = wobble)",
                  fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Substitution direction", fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


def plot_substitution_summary(sub_summary, out_path: Path, dpi: int = 300) -> Path:
    """Grouped bar chart: substitution direction × count, colored by synonymous/non-synonymous.

    Replaces the two-panel heatmap approach with a single figure that directly
    compares synonymous vs non-synonymous counts per substitution direction.
    Counts are pooled across codon positions (1–3) and sorted by total count.
    """
    apply_theme()
    if sub_summary.empty:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "No substitutions observed.", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)

    df = sub_summary.copy()
    df["Direction"] = df["Yeast_Base"] + "→" + df["Human_Base"]
    df["Type"] = df["Same_AA"].map({True: "Synonymous", False: "Non-synonymous"})

    # Pool counts across codon positions; sort directions by total descending
    agg = df.groupby(["Direction", "Type"])["Count"].sum().reset_index()
    total_order = agg.groupby("Direction")["Count"].sum().sort_values(ascending=False).index
    directions = list(total_order)

    syn_map = (agg[agg["Type"] == "Synonymous"]
               .set_index("Direction")["Count"]
               .reindex(directions, fill_value=0))
    nonsyn_map = (agg[agg["Type"] == "Non-synonymous"]
                  .set_index("Direction")["Count"]
                  .reindex(directions, fill_value=0))

    x = np.arange(len(directions))
    width = 0.38

    fig, ax = plt.subplots(figsize=(11, 5))

    bars_syn = ax.bar(x - width / 2, syn_map.values, width,
                      label="Synonymous", color="#2166AC", alpha=0.85)
    bars_non = ax.bar(x + width / 2, nonsyn_map.values, width,
                      label="Non-synonymous", color="#D6604D", alpha=0.85)

    # Value labels above each bar
    for bar in list(bars_syn) + list(bars_non):
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.15,
                    str(int(h)), ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(directions, fontsize=TICK_FONTSIZE)
    ax.set_ylabel("Total substitution count (positions 1–3 pooled)", fontsize=LABEL_FONTSIZE)
    ax.set_xlabel("Substitution direction  (Yeast base → Human base)", fontsize=LABEL_FONTSIZE)
    ax.legend(fontsize=TICK_FONTSIZE, frameon=True, framealpha=0.9)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)
    ax.margins(y=0.12)

    fig.suptitle(
        "Yeast → Human COX1 substitutions: synonymous vs non-synonymous\n"
        "Codon positions 1–3 pooled;  sorted by total count",
        fontsize=TITLE_FONTSIZE - 1, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)
