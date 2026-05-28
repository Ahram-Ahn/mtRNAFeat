"""Foldedness × normalized-MFE landscape and pairing-bias plots.

Publication aesthetics:
- larger figure, larger fonts (set in viz/style.py)
- KDE contours toned down (alpha + level count)
- experimental scatter labelled with force-repelled labels + leader lines
- legend outside the axes so it doesn't sit on top of points
- separate per-species panels for the overlay so Human/Yeast labels don't collide
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from mtrnafeat.viz.samples import ordered_samples, sample_color, sample_palette
from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    LINEWIDTH,
    TITLE_FONTSIZE,
    apply_theme,
    legend_outside,
    panel_label,
    repel_labels,
    style_axis,
)


def landscape_overlay(sim_df, exp_df, out_path: Path, dpi: int = 300) -> Path:
    """Two-panel figure: Human (left) and Yeast (right), each showing the
    species-relevant simulated KDE clouds + the experimental transcripts.

    Each panel is filtered to its species' conditions only (Sim Yeast clouds
    are not shown in the Human panel, and vice versa). The empirical cloud is
    color-filled as a gradient; symmetric-GC reference clouds are drawn as
    distinct contour outlines so the empirical cloud stays the visual focus.
    """
    apply_theme()
    species_list = ordered_samples(exp_df["Species"].unique())
    n_sp = max(len(species_list), 1)
    cols = min(3, n_sp)
    rows = math.ceil(n_sp / cols)
    fig, axes_grid = plt.subplots(rows, cols, figsize=(7.6 * cols, 6.3 * rows), sharey=True, squeeze=False)
    axes = list(axes_grid.flat)
    for ax in axes[n_sp:]:
        ax.axis("off")
    for ax, letter in zip(axes[:n_sp], "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[:n_sp], strict=True):
        panel_label(ax, letter)

    species_palette = sample_palette(species_list)
    # Distinct hues for filled overlapping KDEs. Empirical = blue (primary),
    # symmetric-GC references in warm/contrasting hues so overlaps remain
    # readable. Alpha kept moderate so blended regions still show through.
    empirical_fill = "#3B7FB7"          # saturated blue
    ref_fill_palette = ["#E68A2E", "#7E57C2", "#3CA46B"]  # orange / purple / green

    for ax, species in zip(axes[:n_sp], species_list, strict=True):
        sp_mask = sim_df["Condition"].str.startswith(f"Sim {species}")
        sp_sim = sim_df[sp_mask]
        conditions = list(sp_sim["Condition"].unique())
        # Empirical condition is the one with per-nucleotide-frequency label
        # (e.g. "Sim Human (0.32A 0.25U 0.12G 0.31C)"); the others are
        # symmetric-GC references (e.g. "Sim Human (46% GC)").
        empirical_conds = [c for c in conditions if "GC)" not in c]
        ref_conds = [c for c in conditions if "GC)" in c]

        panel_handles: list = []
        panel_labels: list[str] = []

        for cond in empirical_conds:
            sub = sp_sim[sp_sim["Condition"] == cond]
            try:
                sns.kdeplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, fill=True, alpha=0.50, color=empirical_fill,
                            levels=6, thresh=0.05, warn_singular=False)
            except Exception:
                sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                                ax=ax, color=empirical_fill, alpha=0.50, s=12)
            panel_handles.append(Patch(facecolor=empirical_fill, alpha=0.55))
            panel_labels.append(cond)

        for j, cond in enumerate(ref_conds):
            sub = sp_sim[sp_sim["Condition"] == cond]
            color = ref_fill_palette[j % len(ref_fill_palette)]
            try:
                sns.kdeplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, fill=True, alpha=0.40, color=color,
                            levels=5, thresh=0.08, warn_singular=False)
            except Exception:
                sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                                ax=ax, color=color, alpha=0.40, s=10)
            panel_handles.append(Patch(facecolor=color, alpha=0.45))
            panel_labels.append(cond)

        sub = exp_df[exp_df["Species"] == species]
        if not sub.empty:
            sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, color=species_palette[species],
                            s=140, edgecolor="black", linewidth=1.4, zorder=5)
            ax.margins(x=0.16, y=0.28)
            repel_labels(ax,
                         xs=sub["Normalized_MFE_per_nt"].values,
                         ys=sub["Foldedness_Pct"].values,
                         labels=sub["Gene"].values,
                         color=species_palette[species], fontsize=8, k=42)
            panel_handles.append(Line2D([0], [0], marker="o", color="w",
                                        markerfacecolor=species_palette[species],
                                        markersize=11, markeredgecolor="black"))
            panel_labels.append(f"{species} (in vivo)")

        ax.set_title(f"{species} mt-mRNA", fontsize=TITLE_FONTSIZE, pad=10)
        ax.set_xlabel(r"Normalized MFE  ($\Delta$G kcal/mol per nt)", fontsize=LABEL_FONTSIZE)
        if ax is axes[0]:
            ax.set_ylabel("Structured percentage (%)", fontsize=LABEL_FONTSIZE)
        else:
            ax.set_ylabel("")
        ax.margins(x=0.16, y=0.28)
        style_axis(ax)
        legend_outside(ax, handles=panel_handles, labels=panel_labels,
                       position="bottom", frameon=True, framealpha=0.92,
                       fontsize=8, ncol=1)

    fig.suptitle("In vivo mt-mRNA structures vs. simulated thermodynamic null",
                 fontsize=TITLE_FONTSIZE + 1, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
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
    species_palette = sample_palette([species])

    sim_sub = sim_df[sim_df["Species"] == species]
    sim_conditions = list(sim_sub["Condition"].unique())

    # Reference contours from the symmetric-GC clouds (Species == "n/a")
    refs = _SPECIES_REFERENCE_CONTOURS.get(species, ())
    ref_sub = sim_df[sim_df["Condition"].isin(refs)] if refs else sim_df.iloc[0:0]

    # All KDEs are color-filled (no outline-only contours). Distinct hues +
    # moderate alpha so the empirical cloud and symmetric-GC references stay
    # visually separable even when they overlap.
    empirical_color = "#3B7FB7"  # saturated blue (primary)
    ref_colors = ["#E68A2E", "#7E57C2", "#3CA46B"]  # orange / purple / green

    fig, ax = plt.subplots(figsize=(8, 6))

    legend_handles: list = []
    legend_labels: list[str] = []

    for cond in sim_conditions:
        sub = sim_sub[sim_sub["Condition"] == cond]
        try:
            sns.kdeplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                        ax=ax, fill=True, alpha=0.50, color=empirical_color,
                        levels=6, thresh=0.05, warn_singular=False)
        except Exception:
            sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, color=empirical_color, alpha=0.50, s=12)
        legend_handles.append(Patch(facecolor=empirical_color, alpha=0.55))
        legend_labels.append(cond)

    for i, cond in enumerate(refs):
        sub = ref_sub[ref_sub["Condition"] == cond]
        if sub.empty:
            continue
        color = ref_colors[i % len(ref_colors)]
        try:
            sns.kdeplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                        ax=ax, fill=True, alpha=0.40, color=color,
                        levels=5, thresh=0.08, warn_singular=False)
        except Exception:
            sns.scatterplot(data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                            ax=ax, color=color, alpha=0.40, s=10)
        legend_handles.append(Patch(facecolor=color, alpha=0.45))
        legend_labels.append(cond)

    exp_sub = exp_df[exp_df["Species"] == species]
    if not exp_sub.empty:
        color = species_palette.get(species, "#333333")
        sns.scatterplot(data=exp_sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                        ax=ax, color=color, s=140, edgecolor="black",
                        linewidth=1.4, zorder=5)
        ax.margins(x=0.16, y=0.28)
        repel_labels(ax,
                     xs=exp_sub["Normalized_MFE_per_nt"].values,
                     ys=exp_sub["Foldedness_Pct"].values,
                     labels=exp_sub["Gene"].values,
                     color=color, fontsize=8, k=42)
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
    legend_outside(ax, handles=legend_handles, labels=legend_labels,
                   position="right", fontsize=9, frameon=True, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


def landscape_overlay_regions(sim_df, region_df, out_path: Path, species: str,
                               region_null_df=None, dpi: int = 300) -> Path:
    """Single-sample overlay with yeast UTR/CDS null contours and region points."""
    apply_theme()
    sim_sub = sim_df[sim_df["Species"] == species]
    exp_sub = region_df[region_df["Species"] == species]
    null_sub = (
        region_null_df[region_null_df["Species"] == species].copy()
        if region_null_df is not None and not region_null_df.empty
        else None
    )
    if exp_sub.empty:
        return Path(out_path)

    fig, ax = plt.subplots(figsize=(9.4, 6.9))
    handles: list = []
    labels: list[str] = []

    if null_sub is not None and not null_sub.empty:
        null_sub["Region_Group"] = np.where(null_sub["Region"] == "CDS", "CDS", "UTR")
        contour_defs = [
            ("UTR", "#7570B3", "UTR composition null"),
            ("CDS", "#1B9E77", "CDS composition null"),
        ]
        for group, color, label in contour_defs:
            sub = null_sub[null_sub["Region_Group"] == group]
            if sub.empty:
                continue
            try:
                cmap = sns.light_palette(color, as_cmap=True)
                sns.kdeplot(
                    data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                    ax=ax, fill=True, cmap=cmap, alpha=0.52,
                    levels=8, thresh=0.04, warn_singular=False,
                )
            except Exception:
                sns.scatterplot(
                    data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                    ax=ax, color=color, alpha=0.20, s=14,
                )
            handles.append(Patch(facecolor=color, edgecolor=color, alpha=0.42, label=label))
            labels.append(label)
    elif not sim_sub.empty:
        try:
            sns.kdeplot(
                data=sim_sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                ax=ax, fill=True, alpha=0.45, color="#3B7FB7",
                levels=6, thresh=0.05, warn_singular=False,
            )
        except Exception:
            sns.scatterplot(
                data=sim_sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                ax=ax, color="#3B7FB7", alpha=0.35, s=12,
            )
        handles.append(Patch(facecolor="#3B7FB7", alpha=0.50, label="Simulated sample null"))
        labels.append("Simulated sample null")

    region_colors = {"5'UTR": "#7570B3", "CDS": "#1B9E77", "3'UTR/tail": "#D95F02"}
    region_markers = {"5'UTR": "^", "CDS": "o", "3'UTR/tail": "s"}
    for region in ("5'UTR", "CDS", "3'UTR/tail"):
        sub = exp_sub[exp_sub["Region"] == region]
        if sub.empty:
            continue
        color = region_colors[region]
        marker = region_markers[region]
        sns.scatterplot(
            data=sub, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
            ax=ax, color=color, marker=marker, s=120,
            edgecolor="black", linewidth=1.1, zorder=5,
        )
        repel_labels(
            ax,
            xs=sub["Normalized_MFE_per_nt"].values,
            ys=sub["Foldedness_Pct"].values,
            labels=[f"{g} {region}" for g in sub["Gene"].values],
            color=color, fontsize=7, k=28,
        )
        handles.append(Line2D([0], [0], marker=marker, color="w",
                              markerfacecolor=color, markeredgecolor="black",
                              markersize=9, label=region))
        labels.append(region)

    ax.set_title(f"{species} — DMS regions vs region-specific nulls",
                 fontsize=TITLE_FONTSIZE, pad=10)
    ax.set_xlabel(r"Normalized de novo DMS ΔG  (kcal/mol per nt)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Structured percentage (%)", fontsize=LABEL_FONTSIZE)
    ax.margins(x=0.12, y=0.18)
    style_axis(ax)
    legend_outside(ax, handles=handles, labels=labels,
                   position="right", fontsize=9, frameon=True, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


def _add_composition_ratios(df):
    out = df.copy()
    gc = out.get("Pct_G", 0) + out.get("Pct_C", 0)
    au = out.get("Pct_A", 0) + out.get("Pct_U", 0)
    out["C_share_of_GC_Pct"] = np.where(gc > 0, 100.0 * out.get("Pct_C", 0) / gc, np.nan)
    out["A_share_of_AU_Pct"] = np.where(au > 0, 100.0 * out.get("Pct_A", 0) / au, np.nan)
    return out


def yeast_region_mode_panels(region_null_df, region_df, out_path: Path,
                             species: str = "Yeast", dpi: int = 300) -> Path:
    """Yeast UTR/CDS folding modes with region-matched composition nulls.

    Top row: thermodynamic landscape for each region.
    Bottom row: nucleotide-bias contour for the same region, using GC% vs
    C/(G+C). Experimental DMS region points are overlaid in both rows.
    """
    apply_theme()
    null_sub = _add_composition_ratios(region_null_df[region_null_df["Species"] == species])
    exp_sub = _add_composition_ratios(region_df[region_df["Species"] == species])
    if exp_sub.empty:
        return Path(out_path)

    regions = [r for r in ("5'UTR", "CDS", "3'UTR/tail")
               if (r in set(exp_sub["Region"])) or (r in set(null_sub["Region"]))]
    if not regions:
        return Path(out_path)

    region_colors = {"5'UTR": "#7570B3", "CDS": "#1B9E77", "3'UTR/tail": "#D95F02"}
    markers = {"5'UTR": "^", "CDS": "o", "3'UTR/tail": "s"}
    fig, axes = plt.subplots(2, len(regions), figsize=(5.1 * len(regions), 8.8), squeeze=False)

    for col, region in enumerate(regions):
        color = region_colors.get(region, "#333333")
        marker = markers.get(region, "o")
        sim = null_sub[null_sub["Region"] == region]
        obs = exp_sub[exp_sub["Region"] == region]

        ax = axes[0, col]
        if not sim.empty:
            try:
                sns.kdeplot(
                    data=sim, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                    ax=ax, fill=True, color=color, alpha=0.35,
                    levels=6, thresh=0.05, warn_singular=False,
                )
                sns.kdeplot(
                    data=sim, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                    ax=ax, color=color, linewidths=1.2,
                    levels=5, thresh=0.08, warn_singular=False,
                )
            except Exception:
                sns.scatterplot(
                    data=sim, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                    ax=ax, color=color, alpha=0.20, s=12,
                )
        if not obs.empty:
            sns.scatterplot(
                data=obs, x="Normalized_MFE_per_nt", y="Foldedness_Pct",
                ax=ax, color=color, marker=marker, s=120,
                edgecolor="black", linewidth=1.1, zorder=5,
            )
            repel_labels(
                ax,
                xs=obs["Normalized_MFE_per_nt"].values,
                ys=obs["Foldedness_Pct"].values,
                labels=obs["Gene"].values,
                color=color, fontsize=7, k=24,
            )
        ax.set_title(region, fontsize=TITLE_FONTSIZE - 1, pad=8)
        ax.set_xlabel(r"Region DMS $\Delta$G / nt")
        ax.set_ylabel("Structured percentage (%)" if col == 0 else "")
        ax.margins(x=0.12, y=0.18)
        style_axis(ax)

        ax = axes[1, col]
        sim_comp = sim.dropna(subset=["C_share_of_GC_Pct"])
        obs_comp = obs.dropna(subset=["C_share_of_GC_Pct"])
        if not sim_comp.empty:
            try:
                sns.kdeplot(
                    data=sim_comp, x="Sequence_GC_Pct", y="C_share_of_GC_Pct",
                    ax=ax, color=color, linewidths=1.4,
                    levels=6, thresh=0.05, warn_singular=False,
                )
            except Exception:
                sns.scatterplot(
                    data=sim_comp, x="Sequence_GC_Pct", y="C_share_of_GC_Pct",
                    ax=ax, color=color, alpha=0.18, s=10,
                )
        if not obs_comp.empty:
            sns.scatterplot(
                data=obs_comp, x="Sequence_GC_Pct", y="C_share_of_GC_Pct",
                ax=ax, color=color, marker=marker, s=105,
                edgecolor="black", linewidth=1.0, zorder=5,
            )
        ax.axhline(50, color="#777777", linestyle=":", linewidth=1.0, alpha=0.7)
        ax.set_xlabel("Region GC content (%)")
        ax.set_ylabel("C / (G + C) (%)" if col == 0 else "")
        ax.set_ylim(-3, 103)
        ax.margins(x=0.10, y=0.08)
        style_axis(ax)

    handles = [
        Patch(facecolor=region_colors[r], alpha=0.35, label=f"{r} null contour")
        for r in regions
    ]
    handles.extend([
        Line2D([0], [0], marker=markers.get(r, "o"), color="w",
               markerfacecolor=region_colors[r], markeredgecolor="black",
               markersize=9, label=f"{r} DMS region")
        for r in regions
    ])
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01),
               ncol=min(3, len(handles)), frameon=False, fontsize=9)
    fig.suptitle(
        f"{species} region-specific folding modes\n"
        "Null contours preserve each region's empirical nucleotide composition; dots are genes",
        fontsize=TITLE_FONTSIZE, y=1.01,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 1])
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
    exp_color = sample_color(species)

    bg = biased_gradient_df[biased_gradient_df.get("Species", "") == species]
    exp_sub = exp_df[exp_df["Species"] == species]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    axes = axes.flatten()

    for ax, base in zip(axes, ["A", "C", "G", "U"], strict=True):
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
            ax.margins(x=0.06, y=0.18)
            repel_labels(
                ax,
                xs=exp_sub["Sequence_GC_Pct"].values,
                ys=exp_sub[y_col].values,
                labels=exp_sub["Gene"].values,
                color=exp_color, fontsize=8, k=30,
            )

        ax.set_title(f"{_BASE_LABEL[base]} composition",
                      fontsize=TITLE_FONTSIZE - 1, fontweight="bold", pad=6)
        ax.set_xlabel("Linear sequence GC content (%)", fontsize=LABEL_FONTSIZE)
        ax.set_ylabel(f"{_BASE_LABEL[base]} content (%)", fontsize=LABEL_FONTSIZE)
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.45)
        ax.set_axisbelow(True)
        ax.margins(x=0.06, y=0.18)
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


def _add_paired_nt_columns(df):
    """Add Paired_*_nt_Pct columns (% of sequence in each pair type).

    These equal ``Foldedness_Pct * Paired_*_Pct / 100`` because the existing
    ``Paired_*_Pct`` fields are pair-type ratios (e.g. G-C pairs / total pairs).
    Multiplying by foldedness converts the ratio back into a sequence fraction
    — this is what captures the foldedness drop when G:C is imbalanced.
    """
    out = df.copy()
    if "Foldedness_Pct" in out.columns:
        for pt in ("GC", "AU", "GU"):
            ratio_col = f"Paired_{pt}_Pct"
            if ratio_col in out.columns:
                out[f"Paired_{pt}_nt_Pct"] = out["Foldedness_Pct"] * out[ratio_col] / 100.0
    return out


def paired_nt_fractions(symmetric_gradient_df, biased_gradient_df, exp_df,
                         out_path: Path, species: str,
                         species_freqs: dict | None = None,
                         dpi: int = 300) -> Path:
    """Three-panel figure showing how H-strand-biased composition shifts the
    *absolute* paired-nucleotide budget (not just the pair-type ratio).

    Panel A: overall foldedness (% nt paired) — symmetric vs H-strand.
    Panel B: G-C paired nt as % of sequence + theoretical ceiling 2·min(G,C).
    Panel C: A-U paired nt as % of sequence + theoretical ceiling 2·min(A,U).

    The ceilings make the mechanism explicit: when one of {G,C} is the limiting
    nucleotide, only ``2·min(G,C)/L`` of the sequence can ever sit in G-C
    pairs, so foldedness must drop relative to a symmetric composition at the
    same overall GC%.
    """
    apply_theme()
    exp_color = sample_color(species)

    sym_df = _add_paired_nt_columns(symmetric_gradient_df)
    bias_df = _add_paired_nt_columns(biased_gradient_df)
    bias_sub = bias_df[bias_df.get("Species", "") == species]
    exp_full = _add_paired_nt_columns(exp_df)
    exp_sub = exp_full[exp_full["Species"] == species]

    if species_freqs:
        gc_total = species_freqs["G"] + species_freqs["C"]
        au_total = species_freqs["A"] + species_freqs["U"]
        c_share = species_freqs["C"] / gc_total if gc_total > 0 else 0.5
        a_share = species_freqs["A"] / au_total if au_total > 0 else 0.5
        # Ceilings (% of sequence). With overall GC=x∈[0,1]:
        #   G fraction = (1-c_share)*x, C fraction = c_share*x.
        #   Max GC-pair nt fraction = 2 * min(G,C) = 2 * min(c_share, 1-c_share) * x.
        # AU analogous with a_share against (1 - overall GC).
        gc_ceiling_slope = 2.0 * min(c_share, 1.0 - c_share)  # × GC%
        au_ceiling_slope_au_pct = 2.0 * min(a_share, 1.0 - a_share)  # × (100 - GC)%
        subtitle = (f"H-strand bias: C/(G+C) = {c_share:.2f}, "
                    f"A/(A+U) = {a_share:.2f}")
    else:
        gc_ceiling_slope = 1.0  # symmetric fallback
        au_ceiling_slope_au_pct = 1.0
        subtitle = "Symmetric composition assumed"

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    panel_defs = [
        ("Foldedness_Pct", "Foldedness (% nt paired)", "A. Overall foldedness"),
        ("Paired_GC_nt_Pct", "G-C paired nt (% of sequence)", "B. G-C paired nucleotides"),
        ("Paired_AU_nt_Pct", "A-U paired nt (% of sequence)", "C. A-U paired nucleotides"),
    ]

    for ax, (y_col, ylabel, title) in zip(axes, panel_defs, strict=True):
        sns.lineplot(
            data=sym_df, x="Sequence_GC_Pct", y=y_col, ax=ax,
            color="#b0b0b0", linestyle="--", errorbar=None,
            linewidth=LINEWIDTH * 0.85,
        )
        sns.lineplot(
            data=bias_sub, x="Sequence_GC_Pct", y=y_col, ax=ax,
            color="#333333", errorbar="sd", linewidth=LINEWIDTH,
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
        # Ceiling overlays for panels B and C.
        if y_col == "Paired_GC_nt_Pct":
            xs = [0.0, 100.0]
            sym_ceiling = [0.0, 100.0]                    # symmetric: 2·min = GC%
            bias_ceiling = [0.0, gc_ceiling_slope * 100]  # H-strand: 2·min(c,1-c)·GC%
            ax.plot(xs, sym_ceiling, color="#888888", linestyle=":",
                    linewidth=1.2, alpha=0.6)
            ax.plot(xs, bias_ceiling, color="#000000", linestyle=":",
                    linewidth=1.4, alpha=0.7)
        elif y_col == "Paired_AU_nt_Pct":
            xs = [0.0, 100.0]
            sym_ceiling = [100.0, 0.0]
            bias_ceiling = [au_ceiling_slope_au_pct * 100, 0.0]
            ax.plot(xs, sym_ceiling, color="#888888", linestyle=":",
                    linewidth=1.2, alpha=0.6)
            ax.plot(xs, bias_ceiling, color="#000000", linestyle=":",
                    linewidth=1.4, alpha=0.7)

        ax.set_title(title, fontsize=TITLE_FONTSIZE - 1, pad=8)
        ax.set_xlabel("Linear sequence GC content (%)")
        ax.set_ylabel(ylabel)
        ax.margins(x=0.04, y=0.10)
        style_axis(ax)
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    # Shared legend on the leftmost panel (where data leaves room at upper-left
    # for foldedness which falls/rises modestly).
    legend_handles = [
        Line2D([0], [0], color="#b0b0b0", linestyle="--",
               linewidth=LINEWIDTH * 0.85, label="Symmetric baseline (G=C, A=U)"),
        Line2D([0], [0], color="#333333", linewidth=LINEWIDTH,
               label=f"H-strand baseline ({species})"),
        Line2D([0], [0], color="#888888", linestyle=":", linewidth=1.2,
               label="Ceiling: symmetric"),
        Line2D([0], [0], color="#000000", linestyle=":", linewidth=1.4,
               label="Ceiling: H-strand (2·min)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=exp_color,
               markersize=10, markeredgecolor="black", label=f"Exp: {species}"),
    ]
    legend_outside(axes[-1], handles=legend_handles, labels=[h.get_label() for h in legend_handles],
                   position="right", frameon=True, framealpha=0.92, fontsize=9)

    fig.suptitle(
        f"{species} — H-strand bias reduces foldedness via composition ceilings\n{subtitle}",
        fontsize=TITLE_FONTSIZE, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


def pairing_bias_species(biased_gradient_df, exp_df, out_path: Path, species: str,
                         y_col: str, ylabel: str, include_yx_line: bool = True,
                         symmetric_gradient_df=None, species_freqs: dict | None = None,
                         dpi: int = 300) -> Path:
    """Single-species pairing-bias plot with an H-strand-biased baseline.

    The dark baseline preserves the species-specific C/(G+C) and A/(A+U)
    ratios (computed from H-strand transcripts) as overall GC% sweeps 0→100.
    A lighter dashed curve overlays the symmetric (G=C, A=U) baseline from
    ``symmetric_gradient_df`` if provided, so the H-strand shift is directly
    visible. Title carries the actual ratio values used.
    """
    apply_theme()
    exp_color = sample_color(species)

    bg = biased_gradient_df[biased_gradient_df.get("Species", "") == species]
    exp_sub = exp_df[exp_df["Species"] == species]

    fig, ax = plt.subplots(figsize=(10, 7))
    handles: list = []

    if symmetric_gradient_df is not None and y_col in symmetric_gradient_df.columns:
        sns.lineplot(
            data=symmetric_gradient_df, x="Sequence_GC_Pct", y=y_col, ax=ax,
            color="#b0b0b0", linestyle="--", errorbar=None, linewidth=LINEWIDTH * 0.8,
        )
        handles.append(Line2D([0], [0], color="#b0b0b0", linestyle="--",
                                linewidth=LINEWIDTH * 0.8,
                                label="Symmetric baseline (G=C, A=U)"))

    if not bg.empty and y_col in bg.columns:
        sns.lineplot(
            data=bg, x="Sequence_GC_Pct", y=y_col, ax=ax,
            color="#333333", errorbar="sd", linewidth=LINEWIDTH,
        )
    handles.append(Line2D([0], [0], color="#333333", linewidth=LINEWIDTH,
                            label=f"H-strand baseline ({species})"))

    if include_yx_line:
        ax.plot([0, 100], [0, 100], color="black", linestyle=":", alpha=0.35)
        handles.append(Line2D([0], [0], color="black", linestyle=":", alpha=0.6,
                                label="y = x"))

    if not exp_sub.empty and y_col in exp_sub.columns:
        sns.scatterplot(data=exp_sub, x="Sequence_GC_Pct", y=y_col, ax=ax,
                        color=exp_color, s=150, edgecolor="black",
                        linewidth=1.4, zorder=5)
        repel_labels(ax,
                     xs=exp_sub["Sequence_GC_Pct"].values,
                     ys=exp_sub[y_col].values,
                     labels=exp_sub["Gene"].values,
                     color=exp_color, fontsize=10)
        handles.append(Line2D([0], [0], marker="o", color="w",
                                markerfacecolor=exp_color, markersize=11,
                                markeredgecolor="black",
                                label=f"Exp: {species}"))

    if species_freqs:
        gc = species_freqs["G"] + species_freqs["C"]
        au = species_freqs["A"] + species_freqs["U"]
        c_share = species_freqs["C"] / gc if gc > 0 else 0.5
        a_share = species_freqs["A"] / au if au > 0 else 0.5
        subtitle = (f"H-strand bias: C/(G+C) = {c_share:.2f}, "
                    f"A/(A+U) = {a_share:.2f}")
    else:
        subtitle = "H-strand bias preserved across GC sweep"

    ax.set_title(f"{species} — H-strand-biased baseline vs. experimental\n{subtitle}",
                 fontsize=TITLE_FONTSIZE, pad=10)
    ax.set_xlabel("Linear sequence GC content (%)")
    ax.set_ylabel(ylabel)
    ax.margins(x=0.06, y=0.10)
    style_axis(ax)
    leg = ax.get_legend()
    if leg is not None:
        leg.remove()
    # GC baseline rises from origin (upper-left empty); AU baseline falls from
    # top-left (lower-left empty). Pick a fixed empty corner per pair type to
    # avoid the legend overflowing the axes.
    legend_loc = "lower left" if "A-U" in ylabel else "upper left"
    ax.legend(handles=handles, loc=legend_loc, frameon=True, framealpha=0.92,
              fontsize=10, borderaxespad=0.6)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


def pairing_bias(gradient_df, exp_df, out_path: Path, y_col: str, ylabel: str,
                 include_yx_line: bool = True, dpi: int = 300) -> Path:
    """Combined-species pairing-bias plot with the SYMMETRIC baseline (G=C, A=U).

    Pair with :func:`pairing_bias_species` for the per-species H-strand-biased
    counterpart. Both species' experimental transcripts overlay the single
    symmetric baseline.
    """
    apply_theme()
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.lineplot(data=gradient_df, x="Sequence_GC_Pct", y=y_col, ax=ax,
                 color="#333333", errorbar="sd", linewidth=LINEWIDTH)
    if include_yx_line:
        ax.plot([0, 100], [0, 100], color="black", linestyle=":", alpha=0.35)
    palette = sample_palette(exp_df["Species"].unique())
    handles = [Line2D([0], [0], color="#333333", linewidth=LINEWIDTH,
                       label="Symmetric baseline (G=C, A=U)")]
    if include_yx_line:
        handles.append(Line2D([0], [0], color="black", linestyle=":", alpha=0.6,
                                label="y = x"))
    for sp in ordered_samples(exp_df["Species"].unique()):
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
    ax.set_title("Human + Yeast — symmetric baseline (G=C, A=U)",
                 fontsize=TITLE_FONTSIZE, pad=10)
    ax.set_xlabel("Linear sequence GC content (%)")
    ax.set_ylabel(ylabel)
    ax.margins(x=0.06, y=0.10)
    style_axis(ax)
    legend_loc = "lower left" if "A-U" in ylabel else "upper left"
    ax.legend(handles=handles, loc=legend_loc, frameon=True, framealpha=0.92,
              fontsize=10, borderaxespad=0.6)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)
