"""CoFold parameter-sweep figures: per-gene heatmaps + per-window correlation curves.

One figure per species (Human / Yeast). Mixing species in a single grid was
visually noisy because each species's mt-genes have very different gene-length
distributions and thus very different |CoFold − DMS| magnitudes.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from mtrnafeat.viz.style import apply_theme, style_axis, LABEL_FONTSIZE, TITLE_FONTSIZE


_SPECIES_ORDER = ["Human", "Yeast"]


def _gap_grid_for_species(species_full: pd.DataFrame, species: str,
                            out_path: Path, dpi: int) -> Path:
    apply_theme()
    if species_full.empty:
        return out_path
    genes = sorted(species_full["Gene"].unique())
    n = len(genes)
    cols = min(4, n)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4.0 * cols, 3.4 * rows), squeeze=False)
    for ax in axes.flat[n:]:
        ax.axis("off")
    for ax, gene in zip(axes.flat, genes):
        sub = species_full[species_full["Gene"] == gene].copy()
        pivot = sub.pivot_table(index="alpha", columns="tau", values="Abs_Gap")
        pivot = pivot.sort_index(ascending=False)
        sns.heatmap(pivot, ax=ax, cmap="viridis_r", annot=True, fmt=".2f",
                    cbar_kws={"label": "|CoFold − DMS| (kcal/mol)"}, linewidths=0.4)
        valid = sub.dropna(subset=["Abs_Gap"])
        if not valid.empty:
            best = valid.sort_values("Abs_Gap").iloc[0]
            try:
                row_idx = pivot.index.get_loc(best["alpha"])
                col_idx = pivot.columns.get_loc(best["tau"])
                ax.add_patch(plt.Rectangle((col_idx, row_idx), 1, 1, fill=False,
                                             edgecolor="red", lw=2.0))
            except Exception:
                pass
        ax.set_title(gene, fontsize=TITLE_FONTSIZE - 3)
        ax.set_xlabel("τ (decay nt)", fontsize=LABEL_FONTSIZE - 2)
        ax.set_ylabel("α (kcal/mol)", fontsize=LABEL_FONTSIZE - 2)
    fig.suptitle(f"{species} — CoFold parameter sweep, gap to DMS-eval ΔG (red = best)",
                 fontsize=TITLE_FONTSIZE - 1, y=1.005)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


def gap_heatmap_panels(full: pd.DataFrame, out_dir: Path, plot_format: str,
                        dpi: int = 300) -> list[Path]:
    """Per-gene |CoFold − DMS| heatmap, one figure per species."""
    if full.empty:
        return []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = plot_format.lstrip(".")
    species_present = [s for s in _SPECIES_ORDER if s in full["Species"].unique()]
    if not species_present:
        species_present = sorted(full["Species"].unique())
    paths = []
    for sp in species_present:
        out_path = out_dir / f"cofold_gap_heatmap_{sp.lower()}.{fmt}"
        _gap_grid_for_species(full[full["Species"] == sp], sp, out_path, dpi)
        paths.append(out_path)
    return paths


def _corr_grid_for_species(species_win: pd.DataFrame, species: str,
                             out_path: Path, dpi: int) -> Path:
    apply_theme()
    if species_win.empty:
        return out_path
    genes = sorted(species_win["Gene"].unique())
    n = len(genes)
    cols = min(4, n)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4.4 * cols, 3.4 * rows), squeeze=False)
    for ax in axes.flat[n:]:
        ax.axis("off")
    taus = sorted(species_win["tau"].unique())
    cmap = sns.color_palette("plasma", len(taus))
    for ax, gene in zip(axes.flat, genes):
        sub = species_win[species_win["Gene"] == gene]
        for color, tau in zip(cmap, taus):
            line = sub[sub["tau"] == tau].sort_values("alpha")
            ax.plot(line["alpha"], line["Pearson_r_CoFold_vs_DMS"],
                     "-o", color=color, lw=1.5, ms=5, label=f"τ={int(tau)}")
        ax.axhline(0, color="gray", ls="--", lw=0.7)
        ax.set_xlabel("α (kcal/mol)", fontsize=LABEL_FONTSIZE - 2)
        ax.set_ylabel("Pearson r (CoFold vs DMS)", fontsize=LABEL_FONTSIZE - 2)
        ax.set_title(gene, fontsize=TITLE_FONTSIZE - 3)
        style_axis(ax)
        ax.legend(fontsize=8, loc="best", ncol=2)
    fig.suptitle(f"{species} — CoFold parameter sweep, per-window correlation with DMS ΔG",
                 fontsize=TITLE_FONTSIZE - 1, y=1.005)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


def per_window_corr_curves(win: pd.DataFrame, out_dir: Path, plot_format: str,
                             dpi: int = 300) -> list[Path]:
    """Per-gene Pearson r curves (alpha on x, r on y), one figure per species."""
    if win.empty:
        return []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = plot_format.lstrip(".")
    species_present = [s for s in _SPECIES_ORDER if s in win["Species"].unique()]
    if not species_present:
        species_present = sorted(win["Species"].unique())
    paths = []
    for sp in species_present:
        out_path = out_dir / f"cofold_per_window_corr_{sp.lower()}.{fmt}"
        _corr_grid_for_species(win[win["Species"] == sp], sp, out_path, dpi)
        paths.append(out_path)
    return paths
