"""CoFold parameter-sweep figures.

Two outputs per species (Human / Yeast):

1. ``cofold_gap_strip_{species}`` — strip plot of |CoFold − DMS| with one
   dot per (α, τ) sweep cell, grouped by gene. Replaces the prior
   per-gene heatmap grid, which was visually busy and made it hard to
   compare genes against each other. The strip plot trades the (α, τ)
   surface detail (still in the CSV: ``cofold_grid.csv``) for a compact
   cross-gene summary that's easier to read in a paper.
2. ``cofold_per_window_corr_{species}`` — per-gene Pearson r curves
   between CoFold ΔG and DMS-eval ΔG across the α grid, one line per τ.
   Unchanged from the prior layout; the heatmap was the noisy figure,
   not this one.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    TITLE_FONTSIZE,
    apply_theme,
    legend_outside,
    style_axis,
)

_SPECIES_ORDER = ["Human", "Yeast"]


def _gap_strip_for_species(species_full: pd.DataFrame, species: str,
                           out_path: Path, dpi: int) -> Path:
    """Strip plot: x=gene (sorted by best gap), y=|CoFold − DMS|, hue=τ.

    The best-fit dot per gene gets a red ring so the eye lands on it
    without needing a heatmap. The default CoFold parameters
    (α = 0.5, τ = 640) get a black square outline so the user can see
    how far from optimal the published defaults sit on each gene.
    """
    apply_theme()
    if species_full.empty:
        return out_path
    df = species_full.dropna(subset=["Abs_Gap"]).copy()
    if df.empty:
        return out_path
    gene_order = (df.groupby("Gene")["Abs_Gap"].min()
                    .sort_values().index.tolist())
    df["Gene"] = pd.Categorical(df["Gene"], categories=gene_order, ordered=True)
    df = df.sort_values(["Gene", "tau", "alpha"])

    width = max(7.5, 0.85 * len(gene_order) + 3.5)
    fig, ax = plt.subplots(figsize=(width, 5.6))

    sns.stripplot(
        data=df, x="Gene", y="Abs_Gap", hue="tau",
        palette="cividis", size=6.5, alpha=0.85, jitter=0.18,
        dodge=False, ax=ax,
    )

    best_per_gene = df.loc[df.groupby("Gene", observed=True)["Abs_Gap"].idxmin()]
    ax.scatter(
        best_per_gene["Gene"].astype(str),
        best_per_gene["Abs_Gap"].to_numpy(),
        s=210, facecolor="none", edgecolor="#D62728", linewidth=2.0,
        zorder=10, label="best (α, τ)",
    )

    default_mask = (df["alpha"].round(3) == 0.5) & (df["tau"].round(0) == 640)
    if default_mask.any():
        defaults = df.loc[default_mask]
        ax.scatter(
            defaults["Gene"].astype(str),
            defaults["Abs_Gap"].to_numpy(),
            s=120, facecolor="none", edgecolor="#222222",
            marker="s", linewidth=1.4, zorder=9,
            label="default (α=0.5, τ=640)",
        )

    ax.set_xlabel("Gene  (sorted by best |gap|)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("|CoFold − DMS|  ΔG gap (kcal/mol)", fontsize=LABEL_FONTSIZE)
    ax.set_title(
        f"{species} — CoFold parameter-sweep gap to DMS-eval ΔG",
        fontsize=TITLE_FONTSIZE - 1, pad=10, fontweight="bold",
    )
    plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
    ax.grid(True, axis="y", linestyle=":", linewidth=0.6, alpha=0.4)
    ax.set_axisbelow(True)
    legend_outside(
        ax, position="right", fontsize=9, frameon=False,
        title="τ — decay (nt) /\noverlay markers",
    )
    style_axis(ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def gap_strip_panels(full: pd.DataFrame, out_dir: Path, plot_format: str,
                     dpi: int = 300) -> list[Path]:
    """Per-species strip plot of |CoFold − DMS| across the (α, τ) sweep.

    Replaces the previous per-gene heatmap grid (``gap_heatmap_panels``).
    """
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
        out_path = out_dir / f"cofold_gap_strip_{sp.lower()}.{fmt}"
        _gap_strip_for_species(full[full["Species"] == sp], sp, out_path, dpi)
        paths.append(out_path)
    return paths


# Backwards-compatible alias so callers (cofold command, run-all
# pipeline) don't break during the rename.
gap_heatmap_panels = gap_strip_panels


def _corr_grid_for_species(species_win: pd.DataFrame, species: str,
                             out_path: Path, dpi: int) -> Path:
    apply_theme()
    if species_win.empty:
        return out_path
    genes = sorted(species_win["Gene"].unique())
    n = len(genes)
    cols = min(4, n)
    rows = math.ceil(n / cols)
    # Wider per-panel allowance so the outside τ-legend doesn't squeeze the data axis.
    fig, axes = plt.subplots(rows, cols, figsize=(6.0 * cols, 3.4 * rows), squeeze=False)
    for ax in axes.flat[n:]:
        ax.axis("off")
    taus = sorted(species_win["tau"].unique())
    cmap = sns.color_palette("cividis", len(taus))
    for ax, gene in zip(axes.flat, genes):
        sub = species_win[species_win["Gene"] == gene]
        for color, tau in zip(cmap, taus):
            line = sub[sub["tau"] == tau].sort_values("alpha")
            ax.plot(line["alpha"], line["Pearson_r_CoFold_vs_DMS"],
                     "-o", color=color, lw=1.5, ms=5, label=f"τ={int(tau)}")
        ax.axhline(0, color="gray", ls="--", lw=0.7)
        ax.set_xlabel("α — penalty strength (kcal/mol)", fontsize=LABEL_FONTSIZE - 2)
        ax.set_ylabel("Pearson r (CoFold vs DMS)", fontsize=LABEL_FONTSIZE - 2)
        ax.set_title(gene, fontsize=TITLE_FONTSIZE - 3)
        style_axis(ax)
        legend_outside(ax, position="right", fontsize=8, frameon=False,
                       title="τ — decay (nt)", title_fontsize=8)
    fig.suptitle(f"{species} — CoFold parameter sweep, per-window correlation with DMS ΔG\n"
                 "f(d) = α · (1 − exp(−d / τ));  larger α → stronger long-range penalty,  larger τ → penalty kicks in only at long distance",
                 fontsize=TITLE_FONTSIZE - 2, y=1.01)
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
