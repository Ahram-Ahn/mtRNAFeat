"""Yeast↔human COX1 comparative plots: substitution map + ΔG-difference track."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

from mtrnafeat.viz.style import LABEL_FONTSIZE, TITLE_FONTSIZE, apply_theme, style_axis


def plot_dG_track_species(track_df, species: str, out_path: Path, dpi: int = 300) -> Path:
    """Single-species COX1 local-ΔG track (no cross-species comparison)."""
    apply_theme()
    col = "yeast_dG" if species.lower().startswith("yeast") else "human_dG"
    color = "#FF7F0E" if col == "yeast_dG" else "#D62728"
    df = track_df.dropna(subset=[col])
    fig, ax = plt.subplots(figsize=(15, 4.6))
    ax.plot(df["col"], df[col], color=color, lw=1.6)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.7)
    ax.set_xlabel("Codon column (PAL2NAL alignment)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"{species} COX1 local ΔG (kcal/mol)", color=color, fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{species} COX1 — local ΔG along the alignment",
                 fontsize=TITLE_FONTSIZE - 2, pad=10)
    style_axis(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)


# Backward-compat alias: the old combined-track function now defaults to
# emitting only the human panel (the comparative-difference panel was
# explicitly requested removed).

def plot_dG_track(track_df, out_path: Path, dpi: int = 300) -> Path:
    """Deprecated: emit per-species figures via plot_dG_track_species instead."""
    return plot_dG_track_species(track_df, "Human", out_path, dpi=dpi)


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
