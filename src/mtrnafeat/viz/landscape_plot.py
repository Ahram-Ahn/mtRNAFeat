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


_SPECIES_REFERENCE_CONTOURS = {
    "Human": ("Sim Human (46% GC)",),
    "Yeast": ("Sim Yeast 5' UTR (7% GC)", "Sim Yeast CDS (30% GC)"),
}


def landscape_overlay_one(sim_df, exp_df, out_path: Path, species: str,
                           dpi: int = 300) -> Path:
    """Single-panel figure for one species: empirical sim cloud + GC-reference
    contours + experimental scatter.

    The empirical (species-specific A/U/G/C) cloud sets the primary contour;
    additional GC-only reference contours are overlaid so the reader can see
    how the empirical cloud compares against well-defined GC isolines (e.g.
    7% GC for yeast 5'UTR, 30% GC for yeast CDS, 46% GC for human).
    """
    apply_theme()
    species_palette = {"Human": "#D62728", "Yeast": "#FF7F0E"}

    sim_sub = sim_df[sim_df["Species"] == species]
    sim_conditions = list(sim_sub["Condition"].unique())

    # Reference contours from the symmetric-GC clouds (Species == "n/a")
    refs = _SPECIES_REFERENCE_CONTOURS.get(species, ())
    ref_sub = sim_df[sim_df["Condition"].isin(refs)] if refs else sim_df.iloc[0:0]

    # Muted sequential palette for empirical cloud; distinct accent colors
    # per GC reference so they don't blur into a single mass.
    empirical_color = "#7E7E7E"  # neutral grey for the per-species empirical cloud
    ref_colors = ["#1f77b4", "#2ca02c", "#9467bd"]  # blue / green / purple

    fig, ax = plt.subplots(figsize=(8, 6))

    legend_handles: list = []
    legend_labels: list[str] = []

    for cond in sim_conditions:
        sub = sim_sub[sim_sub["Condition"] == cond]
        try:
            sns.kdeplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                        ax=ax, fill=True, alpha=0.28, color=empirical_color,
                        levels=5, thresh=0.10, warn_singular=False)
        except Exception:
            sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, color=empirical_color, alpha=0.28, s=12)
        legend_handles.append(Patch(facecolor=empirical_color, alpha=0.40))
        legend_labels.append(cond)

    for i, cond in enumerate(refs):
        sub = ref_sub[ref_sub["Condition"] == cond]
        if sub.empty:
            continue
        color = ref_colors[i % len(ref_colors)]
        try:
            sns.kdeplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                        ax=ax, fill=False, alpha=0.85, color=color,
                        levels=4, thresh=0.10, linewidths=1.8,
                        warn_singular=False)
        except Exception:
            sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, color=color, alpha=0.40, s=10)
        legend_handles.append(Line2D([0], [0], color=color, linewidth=2.0))
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


_BASE_BIAS_COLORS = {
    "A": "#2ca02c",
    "C": "#1f77b4",
    "G": "#9467bd",
    "U": "#d62728",
}
_BASE_LABEL = {"A": "A", "C": "C", "G": "G", "U": "U / T"}


def per_base_composition_bias(biased_gradient_df, exp_df, out_path: Path,
                               species: str, dpi: int = 300) -> Path:
    """Per-nucleotide composition vs sequence GC%, in pairing_bias style.

    Four small panels (A / C / G / U). Each panel:
      * grey line: the species-biased simulation gradient — preserves the
        empirical C:G and A:U ratios as GC% sweeps, so within a given GC
        ratio the per-species nucleotide imbalance is built in.
      * scatter: each experimental gene as a labelled point.

    This replaces the previous violin (which collapsed the gradient into
    a single column and threw away the GC-axis signal).
    """
    apply_theme()
    species_palette = {"Human": "#D62728", "Yeast": "#FF7F0E"}
    exp_color = species_palette.get(species, "#333333")

    bg = biased_gradient_df[biased_gradient_df.get("Species", "") == species]
    exp_sub = exp_df[exp_df["Species"] == species]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    axes = axes.flatten()

    for ax, base in zip(axes, ["A", "C", "G", "U"]):
        y_col = f"Pct_{base}"
        line_color = _BASE_BIAS_COLORS[base]

        if not bg.empty and y_col in bg.columns:
            sns.lineplot(
                data=bg, x="Sequence_GC_Pct", y=y_col, ax=ax,
                color=line_color, errorbar="sd", linewidth=LINEWIDTH,
                label=f"Biased baseline (preserves {species} C:G & A:U)",
            )

        if not exp_sub.empty and y_col in exp_sub.columns:
            sns.scatterplot(
                data=exp_sub, x="Sequence_GC_Pct", y=y_col, ax=ax,
                color=exp_color, s=110, edgecolor="black", linewidth=1.2,
                zorder=5,
            )
            repel_labels(
                ax,
                xs=exp_sub["Sequence_GC_Pct"].values,
                ys=exp_sub[y_col].values,
                labels=exp_sub["Gene"].values,
                color=exp_color, fontsize=9,
            )

        ax.set_title(f"{_BASE_LABEL[base]} composition",
                      fontsize=TITLE_FONTSIZE - 1, fontweight="bold", pad=6)
        ax.set_xlabel("Linear sequence GC content (%)", fontsize=LABEL_FONTSIZE)
        ax.set_ylabel(f"{_BASE_LABEL[base]} content (%)", fontsize=LABEL_FONTSIZE)
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.45)
        ax.set_axisbelow(True)
        ax.margins(x=0.04, y=0.10)
        style_axis(ax)
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    fig.suptitle(
        f"{species} — individual nucleotide composition vs sequence GC%\n"
        f"baseline preserves species-specific C/(G+C) and A/(A+U) bias",
        fontsize=TITLE_FONTSIZE, fontweight="bold", y=1.02,
    )
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
