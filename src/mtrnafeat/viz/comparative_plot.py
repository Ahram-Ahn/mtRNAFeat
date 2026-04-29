"""Yeast↔human COX1 comparative plots: substitution heatmap + directional flux."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

from mtrnafeat.viz.style import LABEL_FONTSIZE, TITLE_FONTSIZE, apply_theme


def plot_directional_flux(flux_df, out_path: Path, dpi: int = 300) -> Path:
    """12 (from→to) × 3 (codon position) heatmap of bias score, BH-significant
    cells outlined in red. Colored by Bias_Score (observed - expected fraction)."""
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
    fig, ax = plt.subplots(figsize=(7.5, 8))
    sns.heatmap(pivot, ax=ax, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                linewidths=0.4, cbar_kws={"label": "Bias score (obs − exp)"})
    # Outline BH-significant cells.
    for i, direction in enumerate(pivot.index):
        for j, pos in enumerate(pivot.columns):
            try:
                q = qpivot.loc[direction, pos]
            except KeyError:
                q = 1.0
            if q < 0.05:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                                            edgecolor="black", lw=2.0))
    ax.set_title("Yeast→Human substitution flux (BH-significant outlined, q<0.05)",
                 fontsize=TITLE_FONTSIZE - 2, pad=10)
    ax.set_xlabel("Codon position (1=1st, 2=2nd, 3=wobble)")
    ax.set_ylabel("Substitution direction")
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def plot_substitution_summary(sub_summary, out_path: Path, dpi: int = 300) -> Path:
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
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, syn_flag, label in [(axes[0], True, "Synonymous"), (axes[1], False, "Non-synonymous")]:
        sub = df[df["Same_AA"] == syn_flag]
        if sub.empty:
            ax.set_title(f"{label} (no events)")
            continue
        pivot = sub.pivot_table(index="Direction", columns="Position",
                                 values="Count", aggfunc="sum", fill_value=0)
        sns.heatmap(pivot, ax=ax, cmap="Blues", annot=True, fmt="d", linewidths=0.5,
                    cbar_kws={"label": "Count"})
        ax.set_title(f"{label} substitutions (Yeast→Human, COX1)", fontsize=12)
        ax.set_xlabel("Codon position (1=1st, 2=2nd, 3=wobble)")
        ax.set_ylabel("Substitution direction")
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
