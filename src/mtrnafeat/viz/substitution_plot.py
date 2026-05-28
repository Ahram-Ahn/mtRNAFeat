"""Visualization of the substitution-thermo permutation test.

Per (species, gene) panel: KDE of three null pools, with two vertical
lines for the wild-type:
  - solid black  : Vienna MFE of the WT sequence (apples-to-apples vs pool)
  - dashed red   : Vienna `eval_structure` of the experimental DMS
                   dot-bracket (read from the .db's structure line, NOT
                   the .db header MFE) on the same chunk length

Per-species figures so human / yeast aren't crammed into one grid.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    LEGEND_FONTSIZE,
    TITLE_FONTSIZE,
    apply_theme,
    style_axis,
)

_POOL_COLORS = {
    "flat_acgu": "#D62728",
    "positional_acgu": "#9467BD",
    "synonymous": "#2CA02C",
}
_SPECIES_ORDER = ["Human", "Yeast"]


def _kde_for_species(species_dist: pd.DataFrame, species: str,
                       out_path: Path, dpi: int) -> Path:
    apply_theme()
    genes = sorted(species_dist["Gene"].unique())
    n = len(genes)
    if n == 0:
        fig = plt.figure(figsize=(4, 3))
        plt.text(0.5, 0.5, f"No {species} substitution data", ha="center", va="center")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return out_path
    cols = min(4, n)
    rows = math.ceil(n / cols)
    # One shared legend avoids repeating five entries in every small panel.
    fig, axes = plt.subplots(rows, cols, figsize=(4.9 * cols + 1.9, 3.2 * rows), squeeze=False)
    for ax in axes.flat[n:]:
        ax.axis("off")
    for ax, gene in zip(axes.flat, genes, strict=False):
        sub = species_dist[species_dist["Gene"] == gene]
        wt_mfe_row = sub[sub["Pool"] == "WildType_MFE"]["MFE_kcal_per_mol"]
        wt_dms_row = sub[sub["Pool"] == "WildType_DMS_Eval"]["MFE_kcal_per_mol"]
        wt_mfe = float(wt_mfe_row.iloc[0]) if not wt_mfe_row.empty else None
        wt_dms = float(wt_dms_row.iloc[0]) if not wt_dms_row.empty else float("nan")
        for pool, color in _POOL_COLORS.items():
            vals = sub[sub["Pool"] == pool]["MFE_kcal_per_mol"].values
            if len(vals) > 5:
                sns.kdeplot(vals, ax=ax, color=color, fill=True, alpha=0.3,
                            linewidth=1.6)
        if wt_mfe is not None:
            ax.axvline(wt_mfe, color="black", linestyle="-", linewidth=1.8,
                       label="WT Vienna MFE")
        if wt_dms is not None and np.isfinite(wt_dms):
            ax.axvline(wt_dms, color="#D62728", linestyle="--", linewidth=1.8,
                       label="DMS structure, Vienna ΔG")
        ax.set_title(gene, fontsize=TITLE_FONTSIZE - 3)
        ax.set_xlabel(r"$\Delta$G (kcal/mol)", fontsize=LABEL_FONTSIZE - 2)
        ax.set_ylabel("density", fontsize=LABEL_FONTSIZE - 2)
        style_axis(ax)
        # Ensure both WT vertical lines fall inside the visible x-range.
        # KDE auto-xlim is dominated by the null pools and can be far from
        # the WT (e.g. yeast ATP9: WT_MFE = -68.9 vs pool means around -38).
        x0, x1 = ax.get_xlim()
        wt_xs = [v for v in (wt_mfe, wt_dms) if v is not None and np.isfinite(v)]
        if wt_xs:
            x0 = min(x0, min(wt_xs))
            x1 = max(x1, max(wt_xs))
            pad = 0.05 * (x1 - x0)
            ax.set_xlim(x0 - pad, x1 + pad)
    legend_handles = [
        Patch(facecolor=color, edgecolor=color, alpha=0.30, label=pool)
        for pool, color in _POOL_COLORS.items()
    ]
    legend_handles.extend([
        Line2D([0], [0], color="black", lw=1.8, label="WT Vienna MFE"),
        Line2D([0], [0], color="#D62728", lw=1.8, linestyle="--",
               label="DMS structure, Vienna ΔG"),
    ])
    fig.legend(handles=legend_handles, loc="center right", bbox_to_anchor=(0.985, 0.52),
               frameon=False, fontsize=LEGEND_FONTSIZE + 1)
    fig.suptitle(f"{species} — substitution-thermodynamic permutation test (Vienna MFE)\n"
                 "DMS ΔG = Vienna eval_structure on the .db dot-bracket (header MFE not used)",
                 fontsize=TITLE_FONTSIZE - 2, y=1.01)
    fig.tight_layout(rect=[0, 0, 0.88, 1.0])
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


def kde_panels(dist: pd.DataFrame, out_dir: Path, plot_format: str,
                 dpi: int = 300) -> list[Path]:
    """One KDE panel-grid figure per species."""
    if dist.empty:
        return []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = plot_format.lstrip(".")
    species_present = [s for s in _SPECIES_ORDER if s in dist["Species"].unique()]
    if not species_present:
        species_present = sorted(dist["Species"].unique())
    paths = []
    for sp in species_present:
        out_path = out_dir / f"substitution_kde_panels_{sp.lower()}.{fmt}"
        _kde_for_species(dist[dist["Species"] == sp], sp, out_path, dpi)
        paths.append(out_path)
    return paths


def _zheat_for_species(sub_summary: pd.DataFrame, species: str,
                         out_path: Path, dpi: int) -> Path:
    apply_theme()
    if sub_summary.empty:
        fig = plt.figure(figsize=(4, 3))
        plt.text(0.5, 0.5, f"No {species} data", ha="center", va="center")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return out_path
    pivot = sub_summary.pivot_table(
        index="Gene", columns="Pool", values="Z_WT_MFE_vs_Pool"
    )
    pool_order = ["flat_acgu", "positional_acgu", "synonymous"]
    pivot = pivot.reindex(columns=pool_order)
    # Wider figure to fit the title comfortably; minimum height keeps the
    # heatmap readable even for one-gene smoke runs.
    fig, ax = plt.subplots(figsize=(10.0, max(4.0, 0.45 * len(pivot) + 2.5)))
    sns.heatmap(pivot, ax=ax, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                cbar_kws={"label": "Z(WT_MFE − Pool)", "pad": 0.02}, linewidths=0.4)
    ax.set_title(f"{species} — Wild-type Vienna ΔG vs synonymous-recoding null pools",
                 fontsize=TITLE_FONTSIZE - 2, pad=12)
    ax.set_xlabel("Null pool", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Gene", fontsize=LABEL_FONTSIZE)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


def z_heatmap(summary: pd.DataFrame, out_dir: Path, plot_format: str,
                dpi: int = 300) -> list[Path]:
    """One Z-heatmap per species."""
    if summary.empty:
        return []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = plot_format.lstrip(".")
    species_present = [s for s in _SPECIES_ORDER if s in summary["Species"].unique()]
    if not species_present:
        species_present = sorted(summary["Species"].unique())
    paths = []
    for sp in species_present:
        out_path = out_dir / f"substitution_z_heatmap_{sp.lower()}.{fmt}"
        _zheat_for_species(summary[summary["Species"] == sp], sp, out_path, dpi)
        paths.append(out_path)
    return paths


def _ensure_effect_columns(summary: pd.DataFrame) -> pd.DataFrame:
    out = summary.copy()
    if "Delta_WT_MFE_per_nt_minus_Pool_Mean" not in out.columns:
        out["Delta_WT_MFE_per_nt_minus_Pool_Mean"] = (
            (out["WT_MFE"] - out["Pool_Mean_MFE"]) / out["Length_nt"]
        )
    if "Delta_WT_DMS_per_nt_minus_Pool_Mean" not in out.columns:
        out["Delta_WT_DMS_per_nt_minus_Pool_Mean"] = (
            (out["WT_DMS_Eval"] - out["Pool_Mean_MFE"]) / out["Length_nt"]
        )
    return out


def _effect_shift_for_species(sub_summary: pd.DataFrame, species: str,
                              out_path: Path, dpi: int) -> Path:
    apply_theme()
    sub_summary = _ensure_effect_columns(sub_summary)
    if sub_summary.empty:
        fig = plt.figure(figsize=(4, 3))
        plt.text(0.5, 0.5, f"No {species} data", ha="center", va="center")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return out_path

    pool_order = ["flat_acgu", "positional_acgu", "synonymous"]
    syn = sub_summary[sub_summary["Pool"] == "synonymous"]
    if syn.empty:
        gene_order = sorted(sub_summary["Gene"].unique())
    else:
        gene_order = (
            syn.sort_values("Delta_WT_MFE_per_nt_minus_Pool_Mean")["Gene"]
            .drop_duplicates()
            .tolist()
        )
    y_lookup = {gene: i for i, gene in enumerate(gene_order)}

    fig, axes = plt.subplots(
        1, 3, figsize=(15.8, max(4.5, 0.42 * len(gene_order) + 2.2)),
        sharey=True,
    )
    all_vals = []
    for col in ("Delta_WT_MFE_per_nt_minus_Pool_Mean",
                "Delta_WT_DMS_per_nt_minus_Pool_Mean"):
        vals = sub_summary[col].replace([np.inf, -np.inf], np.nan).dropna().tolist()
        all_vals.extend(vals)
    xmax = max(abs(v) for v in all_vals) if all_vals else 0.02
    xlim = (-1.12 * xmax, 1.12 * xmax)

    for ax, pool in zip(axes, pool_order, strict=True):
        pool_df = sub_summary[sub_summary["Pool"] == pool]
        for _, row in pool_df.iterrows():
            gene = row["Gene"]
            y = y_lookup.get(gene)
            if y is None:
                continue
            wt = float(row["Delta_WT_MFE_per_nt_minus_Pool_Mean"])
            dms = float(row["Delta_WT_DMS_per_nt_minus_Pool_Mean"])
            if np.isfinite(wt) and np.isfinite(dms):
                ax.plot([wt, dms], [y, y], color="#B0B0B0", linewidth=1.0, zorder=1)
            if np.isfinite(wt):
                ax.scatter(wt, y, s=62, color="black", edgecolor="white",
                           linewidth=0.7, zorder=3)
            if np.isfinite(dms):
                ax.scatter(dms, y, s=70, marker="D", color="#D62728",
                           edgecolor="black", linewidth=0.7, zorder=4)
        ax.axvspan(xlim[0], 0, color="#2CA02C", alpha=0.08, zorder=0)
        ax.axvline(0, color="#555555", linestyle=":", linewidth=1.1)
        ax.set_xlim(*xlim)
        ax.set_title(pool, fontsize=TITLE_FONTSIZE - 3, pad=8)
        ax.set_xlabel(r"$\Delta\Delta$G per nt vs pool mean")
        ax.set_yticks(range(len(gene_order)))
        ax.set_ylim(len(gene_order) - 0.5, -0.5)
        if ax is axes[0]:
            ax.set_yticklabels(gene_order)
        else:
            ax.tick_params(labelleft=False)
        style_axis(ax)
    axes[0].set_ylabel("Gene", fontsize=LABEL_FONTSIZE)

    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
               markeredgecolor="white", markersize=8, label="WT Vienna MFE"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#D62728",
               markeredgecolor="black", markersize=8, label="DMS structure ΔG"),
        Patch(facecolor="#2CA02C", alpha=0.08, label="More stable than null mean"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.01),
               ncol=3, frameon=False, fontsize=LEGEND_FONTSIZE)
    fig.suptitle(
        f"{species} — effect size of wild-type sequence and DMS structure vs recoding nulls",
        fontsize=TITLE_FONTSIZE - 2, y=1.02,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def effect_shift_panels(summary: pd.DataFrame, out_dir: Path, plot_format: str,
                        dpi: int = 300) -> list[Path]:
    """Per-species effect-size dot plots.

    Negative ΔΔG means the WT reference is more stable than that null pool's
    mean. This is the reader-facing companion to the KDE panels.
    """
    if summary.empty:
        return []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = plot_format.lstrip(".")
    species_present = [s for s in _SPECIES_ORDER if s in summary["Species"].unique()]
    if not species_present:
        species_present = sorted(summary["Species"].unique())
    paths = []
    for sp in species_present:
        out_path = out_dir / f"substitution_effect_shift_{sp.lower()}.{fmt}"
        _effect_shift_for_species(summary[summary["Species"] == sp], sp, out_path, dpi)
        paths.append(out_path)
    return paths
