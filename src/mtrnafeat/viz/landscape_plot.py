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


def landscape_overlay_one(sim_df, exp_df, out_path: Path, species: str,
                           dpi: int = 300) -> Path:
    """Single-panel figure for one species: simulation KDE cloud + experimental scatter.

    Filters `sim_df` to rows where Species == species (empirical conditions only,
    not symmetric-GC legacy clouds).  Legend is placed inside the axes.
    """
    apply_theme()
    species_palette = {"Human": "#D62728", "Yeast": "#FF7F0E"}

    sim_sub = sim_df[sim_df["Species"] == species]
    sim_conditions = list(sim_sub["Condition"].unique())
    contour_palette = sns.color_palette("viridis", max(len(sim_conditions), 3))

    fig, ax = plt.subplots(figsize=(8, 6))

    legend_handles: list = []
    legend_labels: list[str] = []

    for i, cond in enumerate(sim_conditions):
        sub = sim_sub[sim_sub["Condition"] == cond]
        try:
            sns.kdeplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                        ax=ax, fill=True, alpha=0.30, color=contour_palette[i],
                        levels=4, thresh=0.10, warn_singular=False)
        except Exception:
            sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, color=contour_palette[i], alpha=0.30, s=12)
        legend_handles.append(Patch(facecolor=contour_palette[i], alpha=0.40))
        legend_labels.append(cond)

    exp_sub = exp_df[exp_df["Species"] == species]
    if not exp_sub.empty:
        color = species_palette.get(species, "#333333")
        sns.scatterplot(data=exp_sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                        ax=ax, color=color, s=140, edgecolor="black",
                        linewidth=1.4, zorder=5)
        repel_labels(ax,
                     xs=exp_sub["Normalized_MFE_per_nt"].values,
                     ys=exp_sub["Foldedness_Pct"].values,
                     labels=exp_sub["Gene"].values,
                     color=color, fontsize=11)
        legend_handles.append(Line2D([0], [0], marker="o", color="w",
                                      markerfacecolor=color, markersize=11,
                                      markeredgecolor="black",
                                      label=f"{species} (in vivo)"))
        legend_labels.append(f"{species} (in vivo)")

    ax.set_title(f"{species} mt-mRNA — in vivo vs. simulated null",
                 fontsize=TITLE_FONTSIZE, pad=10)
    ax.set_xlabel(r"Normalized MFE  ($\Delta$G kcal/mol per nt)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Structured percentage (%)", fontsize=LABEL_FONTSIZE)
    ax.margins(x=0.10, y=0.12)
    style_axis(ax)
    ax.legend(handles=legend_handles, labels=legend_labels,
              loc="best", fontsize=10, frameon=True, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


def pairing_bias_species_corrected(species_sim_df, exp_df, out_path: Path,
                                    y_col: str, ylabel: str, species: str,
                                    dpi: int = 300, gene_filter: str | None = None) -> Path:
    """Violin of the species-specific simulation null + experimental gene scatter.

    Uses the empirical per-species nucleotide composition as the reference
    null (passed in as `species_sim_df`), replacing the generic GC-gradient
    baseline.  Experimental genes are shown as a scatter column next to the
    violin.  If `gene_filter` is provided, only that gene's experimental point
    is shown (useful for an ND6-only spotlight panel).
    """
    apply_theme()
    species_palette = {"Human": "#D62728", "Yeast": "#FF7F0E"}
    color = species_palette.get(species, "#333333")

    exp_sub = exp_df[exp_df["Species"] == species].copy()
    if gene_filter is not None:
        exp_sub = exp_sub[exp_sub["Gene"] == gene_filter]

    fig, ax = plt.subplots(figsize=(6, 6))

    sim_vals = species_sim_df[y_col].dropna().values
    if len(sim_vals) > 1:
        parts = ax.violinplot([sim_vals], positions=[0], widths=0.55,
                               showmedians=True, showextrema=True)
        for pc in parts["bodies"]:
            pc.set_facecolor("#7fc7f5")
            pc.set_alpha(0.70)
        parts["cmedians"].set_color("#1565C0")
        parts["cbars"].set_color("#555555")
        parts["cmaxes"].set_color("#555555")
        parts["cmins"].set_color("#555555")

    if not exp_sub.empty:
        y_vals = exp_sub[y_col].values
        jitter = 0.04
        rng = __import__("numpy").random.default_rng(0)
        x_jit = rng.uniform(-jitter, jitter, size=len(y_vals)) + 1
        ax.scatter(x_jit, y_vals, color=color, s=100,
                   edgecolor="black", linewidth=1.2, zorder=5)
        for xi, yi, lbl in zip(x_jit, y_vals, exp_sub["Gene"].values):
            ax.text(xi + 0.06, yi, lbl, fontsize=9, va="center", color="#222222")

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Simulation\n(empirical null)", "Observed\ngenes"],
                        fontsize=LABEL_FONTSIZE - 1)
    ax.set_ylabel(ylabel, fontsize=LABEL_FONTSIZE)
    ax.set_xlim(-0.6, 1.9)
    title = f"{species} — {ylabel}\nvs. nucleotide-corrected null"
    if gene_filter is not None:
        title += f"\n({gene_filter} only)"
    ax.set_title(title, fontsize=TITLE_FONTSIZE - 2, pad=8)
    style_axis(ax)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.45)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
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
