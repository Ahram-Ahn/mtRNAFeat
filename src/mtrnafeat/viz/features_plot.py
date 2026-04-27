"""Heatmap, phase-space contour, base-pair-span boxplot.

Replaces legacy 10's plotting (re-shaped to use the cleaner DataFrames the
new analysis layer produces)."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde

from mtrnafeat.constants import PALETTE
from mtrnafeat.viz.style import apply_theme, style_axis, LABEL_FONTSIZE, TITLE_FONTSIZE


def heatmap_size_ratios(df_motifs: pd.DataFrame, max_size: int, out_path: Path, dpi: int = 300) -> Path:
    """Per-species heatmap of element-size enrichment; Sim columns left, DMS right."""
    apply_theme()
    motifs = ["Macro_Helix", "Hairpin", "Bulge", "Internal_Loop"]
    df = df_motifs[df_motifs["Motif"].isin(motifs)].copy()
    df["Size_Cap"] = np.where(df["Size"] <= max_size, df["Size"], max_size + 1)
    counts = df.groupby(["Species", "Type", "Motif", "Size_Cap"]).size().reset_index(name="Count")
    counts["Total"] = counts.groupby(["Species", "Type", "Motif"])["Count"].transform("sum")
    counts["Ratio_Pct"] = 100.0 * counts["Count"] / counts["Total"]

    species_list = ["Human", "Yeast"]
    fig, axes = plt.subplots(1, len(species_list), figsize=(8 * len(species_list), 7.5))
    if len(species_list) == 1:
        axes = [axes]

    motif_display = [m.replace("_", " ") for m in motifs]
    for ax, species in zip(axes, species_list):
        sub = counts[counts["Species"] == species].copy()
        if sub.empty:
            ax.set_visible(False)
            continue
        sub["Display_Motif"] = sub["Motif"].str.replace("_", " ")
        pivot = sub.pivot_table(index="Size_Cap", columns=["Type", "Display_Motif"],
                                  values="Ratio_Pct", aggfunc="sum").fillna(0)
        ordered_cols = [(t, m) for t in ("Sim", "DMS") for m in motif_display if (t, m) in pivot.columns]
        pivot = pivot.reindex(columns=ordered_cols)
        idx = [str(int(i)) if i != max_size + 1 else f"{max_size}+" for i in pivot.index]
        pivot.index = idx
        col_labels = [m for _, m in ordered_cols]
        pivot.columns = col_labels

        sns.heatmap(pivot, ax=ax, cmap="Reds", annot=False, linewidths=0.5, square=True,
                    cbar_kws={"label": "Enrichment (%)", "shrink": 0.7})
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        ax.tick_params(axis="x", rotation=35, labelsize=11)
        ax.tick_params(axis="y", labelsize=11)
        for lbl in ax.get_xticklabels():
            lbl.set_horizontalalignment("left")
            lbl.set_fontweight("bold")

        sim_count = sum(1 for t, _ in ordered_cols if t == "Sim")
        dms_count = sum(1 for t, _ in ordered_cols if t == "DMS")
        if sim_count:
            ax.annotate("Simulated Baseline", xy=(sim_count / 2.0, -1.4), xycoords="data",
                         ha="center", va="bottom", fontsize=12, weight="bold", color="gray")
            ax.plot([0.1, sim_count - 0.1], [-1.0, -1.0], color="gray", lw=2.5, clip_on=False)
        if dms_count:
            ax.annotate("DMS-Probed (In organello)", xy=(sim_count + dms_count / 2.0, -1.4), xycoords="data",
                         ha="center", va="bottom", fontsize=12, weight="bold", color="black")
            ax.plot([sim_count + 0.1, sim_count + dms_count - 0.1], [-1.0, -1.0],
                     color="black", lw=2.5, clip_on=False)

        ax.set_title(f"{species} Structural Elements: Enrichment Landscape",
                      fontsize=TITLE_FONTSIZE, pad=60, weight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Element size (nt / bp)", fontsize=LABEL_FONTSIZE)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def phase_space(df_motifs: pd.DataFrame, out_path: Path, dpi: int = 300) -> Path:
    """Phase-space scatter with KDE-cloud overlay.

    Uses seaborn's kdeplot like the legacy 10 figure (which renders well on
    this data): `levels=8, alpha=0.5, thresh=0.05`. Falls back to hexbin
    only if seaborn raises (truly degenerate covariance).

    Experimental gene labels are drawn through `repel_labels` (with
    adjustText preferred when available) and have a white-stroke halo so
    they remain legible over the contour fill.
    """
    apply_theme()
    transcripts = []
    for (species, src, gene), grp in df_motifs.groupby(["Species", "Type", "Gene"]):
        avg_stem = grp[grp["Motif"] == "Macro_Helix"]["Size"].mean()
        avg_loop = grp[grp["Motif"].isin(["Hairpin", "Bulge", "Internal_Loop"])]["Size"].mean()
        transcripts.append({
            "Species": species, "Type": src, "Gene": gene,
            "Avg_Macro_Stem": float(avg_stem) if pd.notna(avg_stem) else 0.0,
            "Avg_Total_Loop": float(avg_loop) if pd.notna(avg_loop) else 0.0,
        })
    df_scatter = pd.DataFrame(transcripts)
    species_list = ["Human", "Yeast"]
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), sharex=True, sharey=True)

    if not df_scatter.empty:
        xs_all = df_scatter["Avg_Macro_Stem"].values
        ys_all = df_scatter["Avg_Total_Loop"].values
        xpad = max(0.5, 0.1 * (xs_all.max() - xs_all.min())) if len(xs_all) else 0.5
        ypad = max(0.5, 0.1 * (ys_all.max() - ys_all.min())) if len(ys_all) else 0.5
        x_lim = (xs_all.min() - xpad, xs_all.max() + xpad)
        y_lim = (ys_all.min() - ypad, ys_all.max() + ypad)
    else:
        x_lim, y_lim = (0, 10), (0, 10)

    for ax, sp in zip(axes, species_list):
        sim_sub = df_scatter[(df_scatter["Species"] == sp) & (df_scatter["Type"] == "Sim")]
        exp_sub = df_scatter[(df_scatter["Species"] == sp) & (df_scatter["Type"] == "DMS")]
        cmap = "Reds" if sp == "Human" else "Oranges"
        sim_drawn = False
        if not sim_sub.empty:
            try:
                sns.kdeplot(
                    data=sim_sub, x="Avg_Macro_Stem", y="Avg_Total_Loop",
                    ax=ax, fill=True, cmap=cmap, alpha=0.50,
                    levels=8, thresh=0.05, warn_singular=False,
                )
                sim_drawn = True
            except Exception as e:
                print(f"[features] phase_space: KDE failed for {sp} ({e}); using hexbin.")
                ax.hexbin(sim_sub["Avg_Macro_Stem"], sim_sub["Avg_Total_Loop"],
                           gridsize=18, cmap=cmap, mincnt=1, alpha=0.6)
                sim_drawn = True
        color = PALETTE.get(sp, "red")
        if not exp_sub.empty:
            sns.scatterplot(data=exp_sub, x="Avg_Macro_Stem", y="Avg_Total_Loop",
                            ax=ax, color=color, s=160, edgecolor="black", linewidth=1.4, zorder=5,
                            label="In Organello DMS Transcripts")
            from mtrnafeat.viz.style import repel_labels
            repel_labels(ax,
                          exp_sub["Avg_Macro_Stem"].values, exp_sub["Avg_Total_Loop"].values,
                          exp_sub["Gene"].values, color="black", fontsize=11,
                          use_adjusttext=True)
        if sim_drawn:
            from matplotlib.patches import Patch
            base_color = sns.color_palette(cmap, n_colors=8)[5]
            handles, labels = ax.get_legend_handles_labels()
            handles.insert(0, Patch(facecolor=base_color, alpha=0.5, label="Simulated Landscape"))
            labels.insert(0, "Simulated Landscape")
            ax.legend(handles, labels, loc="upper right", fontsize=11, framealpha=0.9)
        elif not exp_sub.empty:
            ax.legend(loc="upper right", fontsize=11, framealpha=0.9)
        ax.set_xlim(*x_lim)
        ax.set_ylim(*y_lim)
        ax.margins(x=0.10, y=0.12)
        ax.set_title(f"{sp} RNA Stem and Loop Sizes\n(Simulated vs. In Organello DMS)",
                     fontsize=TITLE_FONTSIZE, pad=12)
        ax.set_xlabel("Average Helix Size (nt)", fontsize=LABEL_FONTSIZE)
        ax.set_ylabel("Average Loop Size (nt)", fontsize=LABEL_FONTSIZE)
        style_axis(ax)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)


def span_boxplot(df_spans: pd.DataFrame, out_path: Path, dpi: int = 300) -> Path:
    apply_theme()
    fig, ax = plt.subplots(figsize=(10, 6))
    df = df_spans.copy()
    df["Condition"] = df["Type"].astype(str) + " " + df["Species"].astype(str)
    palette = {"DMS Human": "#D62728", "DMS Yeast": "#FF7F0E",
               "Sim Human": "#F4A6A2", "Sim Yeast": "#FFD7A8"}
    sns.boxplot(data=df, x="Condition", y="Span", ax=ax, hue="Condition",
                palette=palette, width=0.6, fliersize=3, linewidth=1.4, legend=False)
    ax.set_title("Base-pairing span distribution", fontsize=TITLE_FONTSIZE, pad=12)
    ax.set_xlabel("Source / organism", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Pairing distance (nt)", fontsize=LABEL_FONTSIZE)
    style_axis(ax)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
