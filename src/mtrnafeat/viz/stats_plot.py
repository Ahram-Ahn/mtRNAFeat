"""Stats summary visualization: per-species side-by-side boxes."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from mtrnafeat.viz.style import LABEL_FONTSIZE, TITLE_FONTSIZE, apply_theme, style_axis

_SPECIES_PALETTE = {"Human": "#D62728", "Yeast": "#FF7F0E"}


def boxplot_per_species(df_stats: pd.DataFrame, out_path: Path, dpi: int = 300) -> Path:
    """Multi-panel boxplot grid showing whichever per-transcript metrics are
    present in the input frame. Caps at 6 panels so the figure stays readable.

    The current canonical columns are emitted by ``analysis.statistics``:
    Foldedness_Pct, Normalized_MFE_per_nt, Paired_GC_Pct, Paired_AU_Pct,
    Paired_GU_Pct, Sequence_GC_Pct.
    """
    apply_theme()
    metrics: list[tuple[str, str]] = []
    candidates = [
        # Current canonical names
        ("Foldedness_Pct", "Foldedness (%)"),
        ("Normalized_MFE_per_nt", "MFE / nt (kcal/mol)"),
        ("Paired_GC_Pct", "Paired G·C (%)"),
        ("Paired_AU_Pct", "Paired A·U (%)"),
        ("Paired_GU_Pct", "Paired G·U wobble (%)"),
        ("Sequence_GC_Pct", "Sequence GC (%)"),
        # Legacy fallbacks — kept so older CSVs still render
        ("Paired_Fraction", "Paired fraction"),
        ("Mean_BP_Distance", "Mean ensemble distance"),
        ("Mean_BP_Span", "Mean BP span (nt)"),
        ("MFE_per_nt", "MFE / nt (kcal/mol)"),
    ]
    seen: set[str] = set()
    for col, label in candidates:
        if col in df_stats.columns and label not in seen:
            metrics.append((col, label))
            seen.add(label)
    metrics = metrics[:6]
    if not metrics:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No per-transcript metrics available.", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)
    n = len(metrics)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 4.0 * nrows), squeeze=False)
    axes_flat = axes.flatten()
    for ax, (col, label) in zip(axes_flat, metrics):
        sns.boxplot(data=df_stats, x="Species", y=col, ax=ax,
                    hue="Species", palette=_SPECIES_PALETTE, width=0.55,
                    fliersize=3, linewidth=1.4, legend=False)
        sns.stripplot(data=df_stats, x="Species", y=col, ax=ax,
                       color="black", size=3, alpha=0.55, jitter=0.2)
        ax.set_title(label, fontsize=TITLE_FONTSIZE - 3)
        ax.set_xlabel("")
        ax.set_ylabel(label, fontsize=LABEL_FONTSIZE - 1)
        style_axis(ax)
    # Hide any unused axes in the grid
    for ax in axes_flat[len(metrics):]:
        ax.set_visible(False)
    fig.suptitle("Per-transcript structural statistics", fontsize=TITLE_FONTSIZE - 1, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
