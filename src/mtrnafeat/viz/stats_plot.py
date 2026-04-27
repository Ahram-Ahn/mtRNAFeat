"""Stats summary visualization: per-species side-by-side boxes."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from mtrnafeat.viz.style import apply_theme, style_axis, LABEL_FONTSIZE, TITLE_FONTSIZE


_SPECIES_PALETTE = {"Human": "#D62728", "Yeast": "#FF7F0E"}


def boxplot_per_species(df_stats: pd.DataFrame, out_path: Path, dpi: int = 300) -> Path:
    """4-panel boxplot grid: Paired_Fraction, Mean_BP_Distance, Mean_BP_Span, MFE_per_nt."""
    apply_theme()
    metrics = []
    candidates = [
        ("Paired_Fraction", "Paired fraction"),
        ("Mean_BP_Distance", "Mean ensemble distance"),
        ("Mean_BP_Span", "Mean BP span (nt)"),
        ("MFE_per_nt", "MFE / nt (kcal/mol)"),
    ]
    for col, label in candidates:
        if col in df_stats.columns:
            metrics.append((col, label))
    if not metrics:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No per-transcript metrics available.", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4.2))
    if len(metrics) == 1:
        axes = [axes]
    for ax, (col, label) in zip(axes, metrics):
        sns.boxplot(data=df_stats, x="Species", y=col, ax=ax,
                    hue="Species", palette=_SPECIES_PALETTE, width=0.55,
                    fliersize=3, linewidth=1.4, legend=False)
        sns.stripplot(data=df_stats, x="Species", y=col, ax=ax,
                       color="black", size=3, alpha=0.55, jitter=0.2)
        ax.set_title(label, fontsize=TITLE_FONTSIZE - 3)
        ax.set_xlabel("")
        ax.set_ylabel(label, fontsize=LABEL_FONTSIZE - 1)
        style_axis(ax)
    fig.suptitle("Per-transcript structural statistics", fontsize=TITLE_FONTSIZE - 1, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
