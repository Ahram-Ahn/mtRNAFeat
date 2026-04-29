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

from mtrnafeat.viz.style import LABEL_FONTSIZE, TITLE_FONTSIZE, apply_theme, style_axis

_POOL_COLORS = {
    "flat_gc": "#FF7F0E",
    "positional_gc": "#1F77B4",
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
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.2 * rows), squeeze=False)
    for ax in axes.flat[n:]:
        ax.axis("off")
    for ax, gene in zip(axes.flat, genes):
        sub = species_dist[species_dist["Gene"] == gene]
        wt_mfe_row = sub[sub["Pool"] == "WildType_MFE"]["MFE_kcal_per_mol"]
        wt_dms_row = sub[sub["Pool"] == "WildType_DMS_Eval"]["MFE_kcal_per_mol"]
        wt_mfe = float(wt_mfe_row.iloc[0]) if not wt_mfe_row.empty else None
        wt_dms = float(wt_dms_row.iloc[0]) if not wt_dms_row.empty else float("nan")
        for pool, color in _POOL_COLORS.items():
            vals = sub[sub["Pool"] == pool]["MFE_kcal_per_mol"].values
            if len(vals) > 5:
                sns.kdeplot(vals, ax=ax, color=color, fill=True, alpha=0.3,
                            linewidth=1.6, label=pool)
        if wt_mfe is not None:
            ax.axvline(wt_mfe, color="black", linestyle="-", linewidth=1.8,
                       label=f"WT Vienna MFE = {wt_mfe:.1f}")
        if wt_dms is not None and np.isfinite(wt_dms):
            ax.axvline(wt_dms, color="#D62728", linestyle="--", linewidth=1.8,
                       label=f"DMS structure, Vienna ΔG = {wt_dms:.1f}")
        ax.set_title(gene, fontsize=TITLE_FONTSIZE - 3)
        ax.set_xlabel(r"$\Delta$G (kcal/mol)", fontsize=LABEL_FONTSIZE - 2)
        ax.set_ylabel("density", fontsize=LABEL_FONTSIZE - 2)
        style_axis(ax)
        ax.legend(fontsize=8, loc="best")
    fig.suptitle(f"{species} — substitution-thermodynamic permutation test (Vienna MFE)\n"
                 "DMS ΔG = Vienna eval_structure on the .db dot-bracket (header MFE not used)",
                 fontsize=TITLE_FONTSIZE - 2, y=1.01)
    fig.tight_layout()
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
    pivot = pivot.reindex(columns=["flat_gc", "positional_gc", "synonymous"])
    fig, ax = plt.subplots(figsize=(6, max(3.0, 0.35 * len(pivot))))
    sns.heatmap(pivot, ax=ax, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                cbar_kws={"label": "Z(WT_MFE − Pool)"}, linewidths=0.4)
    ax.set_title(f"{species} — Wild-type Vienna ΔG vs synonymous-recoding null pools",
                 fontsize=TITLE_FONTSIZE - 2, pad=10)
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
