"""Foldedness × normalized-MFE landscape and pairing-bias plots.

Publication aesthetics:
- larger figure, larger fonts (set in viz/style.py)
- KDE contours toned down (alpha + level count)
- experimental scatter labelled with force-repelled labels + leader lines
- legend outside the axes so it doesn't sit on top of points
- separate per-species panels for the overlay so Human/Yeast labels don't collide
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    LINEWIDTH,
    TITLE_FONTSIZE,
    apply_theme,
    panel_label,
    repel_labels,
    style_axis,
)


def landscape_overlay(sim_df, exp_df, out_path: Path, dpi: int = 300) -> Path:
    """Two-panel figure: Human (left) and Yeast (right), each showing the
    simulated KDE clouds + the experimental transcripts as scatter."""
    apply_theme()
    species_list = ["Human", "Yeast"]
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
    for ax, letter in zip(axes, ("A", "B")):
        panel_label(ax, letter)

    sim_conditions = list(sim_df["Condition"].unique())
    contour_palette = sns.color_palette("viridis", max(len(sim_conditions), 3))
    species_palette = {"Human": "#D62728", "Yeast": "#FF7F0E"}

    legend_handles: list = []
    legend_labels: list[str] = []

    for ax, species in zip(axes, species_list):
        for i, cond in enumerate(sim_conditions):
            sub = sim_df[sim_df["Condition"] == cond]
            try:
                sns.kdeplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, fill=True, alpha=0.30, color=contour_palette[i],
                            levels=4, thresh=0.10, warn_singular=False)
            except Exception:
                sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                                ax=ax, color=contour_palette[i], alpha=0.30, s=12)
            if ax is axes[0]:
                legend_handles.append(Patch(facecolor=contour_palette[i], alpha=0.40))
                legend_labels.append(cond)

        sub = exp_df[exp_df["Species"] == species]
        if not sub.empty:
            sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, color=species_palette[species],
                            s=140, edgecolor="black", linewidth=1.4, zorder=5)
            repel_labels(ax,
                         xs=sub["Normalized_MFE_per_nt"].values,
                         ys=sub["Foldedness_Pct"].values,
                         labels=sub["Gene"].values,
                         color=species_palette[species], fontsize=11)

        ax.set_title(f"{species} mt-mRNA", fontsize=TITLE_FONTSIZE, pad=10)
        ax.set_xlabel(r"Normalized MFE  ($\Delta$G kcal/mol per nt)", fontsize=LABEL_FONTSIZE)
        if ax is axes[0]:
            ax.set_ylabel("Structured percentage (%)", fontsize=LABEL_FONTSIZE)
        else:
            ax.set_ylabel("")
        ax.margins(x=0.10, y=0.12)
        style_axis(ax)

    legend_handles += [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=species_palette["Human"],
               markersize=11, markeredgecolor="black", label="Human (in vivo)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=species_palette["Yeast"],
               markersize=11, markeredgecolor="black", label="Yeast (in vivo)"),
    ]
    legend_labels += ["Human (in vivo)", "Yeast (in vivo)"]

    fig.legend(legend_handles, legend_labels, loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=True, fontsize=11)
    fig.suptitle("In vivo mt-mRNA structures vs. simulated thermodynamic null",
                 fontsize=TITLE_FONTSIZE + 1, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def gradient_curves(gradient_df, out_path: Path, dpi: int = 300) -> Path:
    apply_theme()
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    sns.lineplot(data=gradient_df, x="GC_Target_Pct", y="Foldedness_Pct",
                 ax=axes[0], color="#1F77B4", errorbar="sd", linewidth=LINEWIDTH)
    axes[0].set_title("A. GC% → Structured percentage", pad=10)
    axes[0].set_xlabel("Sequence GC content (%)")
    axes[0].set_ylabel("Structured percentage (%)")

    sns.lineplot(data=gradient_df, x="GC_Target_Pct", y="Normalized_MFE_per_nt",
                 ax=axes[1], color="#D62728", errorbar="sd", linewidth=LINEWIDTH)
    axes[1].axhline(0, color="black", linestyle="--", alpha=0.4)
    axes[1].set_title("B. GC% → Normalized MFE", pad=10)
    axes[1].set_xlabel("Sequence GC content (%)")
    axes[1].set_ylabel(r"Normalized MFE ($\Delta$G kcal/mol per nt)")

    for ax in axes:
        style_axis(ax)
    fig.suptitle("Continuous thermodynamic gradient (0–100% GC)",
                 fontsize=TITLE_FONTSIZE, y=1.03)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def pairing_bias(gradient_df, exp_df, out_path: Path, y_col: str, ylabel: str,
                 include_yx_line: bool = True, dpi: int = 300) -> Path:
    apply_theme()
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.lineplot(data=gradient_df, x="Sequence_GC_Pct", y=y_col, ax=ax,
                 color="#777777", errorbar="sd", linewidth=LINEWIDTH,
                 label="Thermodynamic baseline (in silico)")
    if include_yx_line:
        ax.plot([0, 100], [0, 100], "k--", alpha=0.4, label="Random assortment (y=x)")
    palette = {"Human": "#D62728", "Yeast": "#FF7F0E"}
    handles = [Line2D([0], [0], color="#777777", linewidth=LINEWIDTH,
                       label="Thermodynamic baseline (in silico)")]
    if include_yx_line:
        handles.append(Line2D([0], [0], color="black", linestyle="--",
                                label="Random assortment (y=x)"))
    for sp in ("Human", "Yeast"):
        sub = exp_df[exp_df["Species"] == sp]
        if sub.empty:
            continue
        sns.scatterplot(data=sub, x="Sequence_GC_Pct", y=y_col, ax=ax,
                        color=palette[sp], s=150, edgecolor="black",
                        linewidth=1.4, zorder=5)
        repel_labels(ax,
                     xs=sub["Sequence_GC_Pct"].values,
                     ys=sub[y_col].values,
                     labels=sub["Gene"].values,
                     color=palette[sp], fontsize=10)
        handles.append(Line2D([0], [0], marker="o", color="w",
                                markerfacecolor=palette[sp], markersize=11,
                                markeredgecolor="black", label=f"Exp: {sp}"))
    ax.set_xlabel("Linear sequence GC content (%)")
    ax.set_ylabel(ylabel)
    ax.margins(x=0.06, y=0.10)
    style_axis(ax)
    ax.legend(handles=handles, bbox_to_anchor=(1.02, 1.0), loc="upper left",
              frameon=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
