"""CoFold parameter-sweep figures.

Outputs:

1. ``cofold_gap_closure`` — Per species (Human | Yeast), plots the
   *fraction of the Vienna→DMS ΔG gap that is closed* as alpha increases,
   with one line per τ and ±1 SD shading across genes. The DMS-target
   reference label sits on the LEFT side so it never sits on top of the
   data lines on the right.

2. ``cofold_parameter_landscape`` — alpha × tau heatmap, color = mean
   |CoFold − DMS| across genes, one panel per species side-by-side. The
   minimum cell is implicit (darkest color); no star marker is drawn so
   the printed values stay legible.

3. Per-gene parameter landscape heatmaps — same alpha × tau grid as the
   summary landscape, but one cell per (species, gene) tile so the reader
   can see which genes drive the species-level optimum.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    TICK_FONTSIZE,
    TITLE_FONTSIZE,
    apply_theme,
    legend_outside,
    panel_label,
    style_axis,
)

_SPECIES_ORDER = ["Human", "Yeast"]
_SPECIES_COLORS = {"Human": "#2166AC", "Yeast": "#D6604D"}


# ---------------------------------------------------------------------------
# Helper: compute fraction-of-gap-closed per (gene, alpha, tau)
# ---------------------------------------------------------------------------

def _add_frac_closed(df: pd.DataFrame) -> pd.DataFrame:
    baseline = (
        df[df["alpha"] == 0.0]
        .groupby(["Species", "Gene"])["CoFold_MFE"]
        .first()
        .rename("Vienna_alpha0")
        .reset_index()
    )
    df = df.merge(baseline, on=["Species", "Gene"], how="left")
    denom = df["DMS_Eval_dG"] - df["Vienna_alpha0"]
    df = df[denom.abs() > 0.01].copy()
    df["frac_closed"] = (df["CoFold_MFE"] - df["Vienna_alpha0"]) / (
        df["DMS_Eval_dG"] - df["Vienna_alpha0"]
    )
    return df


# ---------------------------------------------------------------------------
# Figure 1: gap-closure convergence
# ---------------------------------------------------------------------------

def _draw_closure_panel(ax, species_df: pd.DataFrame, species: str,
                         taus: list, cmap: list) -> None:
    agg = (
        species_df.groupby(["alpha", "tau"])["frac_closed"]
        .agg(["mean", "std"])
        .reset_index()
        .rename(columns={"mean": "mu", "std": "sd"})
    )
    alphas_sorted = sorted(agg["alpha"].unique())

    ax.axhspan(1.0, 1.6, alpha=0.07, color="#CC0000", zorder=0, lw=0)

    for color, tau in zip(cmap, taus):
        sub = agg[agg["tau"] == tau].sort_values("alpha")
        mu = sub["mu"].to_numpy()
        sd = sub["sd"].fillna(0).to_numpy()
        xs = sub["alpha"].to_numpy()
        ax.plot(xs, mu, "-o", color=color, lw=2.0, ms=6,
                label=f"τ = {int(tau)} nt")
        ax.fill_between(xs, mu - sd, mu + sd, alpha=0.13, color=color)

    ax.axhline(1.0, color="#333333", ls="--", lw=1.3, zorder=5)
    # Place the DMS-target text on the LEFT side so it doesn't sit on top
    # of the data lines at high alpha.
    ax.text(
        min(alphas_sorted) + (max(alphas_sorted) - min(alphas_sorted)) * 0.02, 1.035,
        "CoFold = DMS ΔG",
        ha="left", va="bottom", fontsize=TICK_FONTSIZE,
        color="#333333", style="italic",
    )

    ax.set_xlim(min(alphas_sorted) - 0.05, max(alphas_sorted) + 0.05)
    ax.set_ylim(-0.05, 1.55)
    ax.set_xlabel("α — penalty strength (kcal/mol)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(
        "Fraction of ΔG gap closed\n(0 = plain Vienna,  1 = matches DMS)",
        fontsize=LABEL_FONTSIZE,
    )
    ax.set_title(species, fontsize=TITLE_FONTSIZE, fontweight="bold", pad=6)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.55, alpha=0.45)
    ax.set_axisbelow(True)
    style_axis(ax)


def gap_closure_panels(full: pd.DataFrame, out_dir: Path, plot_format: str,
                        dpi: int = 300) -> list[Path]:
    if full.empty:
        return []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = plot_format.lstrip(".")
    out_path = out_dir / f"cofold_gap_closure.{fmt}"

    apply_theme()
    df = _add_frac_closed(full)
    if df.empty:
        return []

    taus = sorted(df["tau"].unique())
    cmap = sns.color_palette("cividis", len(taus))

    species_present = [s for s in _SPECIES_ORDER if s in df["Species"].unique()]
    if not species_present:
        species_present = sorted(df["Species"].unique())

    n_sp = len(species_present)
    fig, axes = plt.subplots(1, n_sp, figsize=(5.5 * n_sp, 4.6), squeeze=False)

    for ax, sp, letter in zip(axes[0], species_present, "AB"):
        _draw_closure_panel(ax, df[df["Species"] == sp], sp, taus, cmap)
        panel_label(ax, letter)

    legend_outside(
        axes[0][-1], position="right", fontsize=TICK_FONTSIZE,
        frameon=False, title="τ — decay length", title_fontsize=TICK_FONTSIZE,
    )

    fig.suptitle(
        "CoFold long-range penalty closes the ΔG gap to DMS-guided structures\n"
        "f(d) = α · (1 − e^{−d/τ});  shading = ±1 SD across genes;  "
        "red zone = overshoot (CoFold more negative than DMS)",
        fontsize=TICK_FONTSIZE, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return [out_path]


# ---------------------------------------------------------------------------
# Figure 2: parameter landscape heatmap (alpha × tau → mean |gap|)
# ---------------------------------------------------------------------------

def _draw_heatmap(ax, pivot: pd.DataFrame, *, vmin: float, vmax: float,
                  title: str, draw_xlabel: bool = True,
                  draw_ylabel: bool = True):
    im = ax.imshow(
        pivot.values,
        aspect="auto",
        cmap="YlOrRd_r",
        vmin=vmin, vmax=vmax,
        origin="upper",
    )
    alphas = [f"{a:.2g}" for a in pivot.columns.tolist()]
    taus = [str(int(t)) for t in pivot.index.tolist()]
    ax.set_xticks(range(len(alphas)))
    ax.set_xticklabels(alphas, fontsize=TICK_FONTSIZE)
    ax.set_yticks(range(len(taus)))
    ax.set_yticklabels(taus, fontsize=TICK_FONTSIZE)
    if draw_xlabel:
        ax.set_xlabel("α — penalty strength (kcal/mol)", fontsize=LABEL_FONTSIZE)
    if draw_ylabel:
        ax.set_ylabel("τ — decay length (nt)", fontsize=LABEL_FONTSIZE)
    ax.set_title(title, fontsize=TITLE_FONTSIZE - 1, fontweight="bold", pad=6)

    # Annotate each cell. Mid-range cells use black; the very darkest cells
    # (small-gap region of the reversed map) use white so the text stays
    # readable on both ends.
    for r in range(pivot.shape[0]):
        for c in range(pivot.shape[1]):
            val = pivot.values[r, c]
            if np.isfinite(val):
                # Reversed YlOrRd: small value -> dark red, large -> pale.
                # Pale background needs dark text; dark background -> light.
                ratio = (val - vmin) / max(1e-9, (vmax - vmin))
                ax.text(c, r, f"{val:.1f}", ha="center", va="center",
                        fontsize=7.0,
                        color="white" if ratio < 0.30 else "#1a1a1a",
                        fontweight="bold")
    return im


def gap_heatmap_panels(full: pd.DataFrame, out_dir: Path, plot_format: str,
                        dpi: int = 300) -> list[Path]:
    """Side-by-side α × τ heatmaps of mean |CoFold − DMS| (Human | Yeast).

    No star marker — the darkest cell is the optimum and the printed
    values are the source of truth.
    """
    if full.empty:
        return []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = plot_format.lstrip(".")
    out_path = out_dir / f"cofold_parameter_landscape.{fmt}"

    apply_theme()

    species_present = [s for s in _SPECIES_ORDER if s in full["Species"].unique()]
    if not species_present:
        species_present = sorted(full["Species"].unique())

    n_sp = len(species_present)
    fig, axes = plt.subplots(1, n_sp, figsize=(4.5 * n_sp, 3.8), squeeze=False)

    pivot_all = (
        full.dropna(subset=["Abs_Gap"])
        .groupby(["Species", "alpha", "tau"])["Abs_Gap"]
        .mean()
        .reset_index()
    )
    vmin = pivot_all["Abs_Gap"].min()
    vmax = pivot_all["Abs_Gap"].max()

    im = None
    for ax, sp, letter in zip(axes[0], species_present, "AB"):
        sub = pivot_all[pivot_all["Species"] == sp]
        pivot = sub.pivot(index="tau", columns="alpha", values="Abs_Gap")
        pivot = pivot.iloc[::-1]
        im = _draw_heatmap(ax, pivot, vmin=vmin, vmax=vmax, title=sp)
        panel_label(ax, letter)

    fig.suptitle(
        "CoFold parameter landscape — mean gap to DMS-evaluated ΔG\n"
        "dark = smaller gap = better fit",
        fontsize=TICK_FONTSIZE, y=1.02,
    )
    fig.tight_layout(rect=[0, 0, 0.88, 1.0])

    if im is not None:
        cbar_ax = fig.add_axes([0.90, 0.15, 0.025, 0.70])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("Mean |CoFold − DMS| ΔG  (kcal/mol)", fontsize=TICK_FONTSIZE)
        cbar.ax.tick_params(labelsize=TICK_FONTSIZE)

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return [out_path]


# Backwards-compatible alias
gap_strip_panels = gap_heatmap_panels


# ---------------------------------------------------------------------------
# Figure 3: per-gene parameter landscape heatmaps
# ---------------------------------------------------------------------------

def per_gene_landscape(full: pd.DataFrame, out_dir: Path, plot_format: str,
                        dpi: int = 300) -> list[Path]:
    """Per-gene α × τ heatmaps of mean |CoFold − DMS|.

    Layout: rows = species, columns = genes within each species, one tile
    per gene with the same colour scale and the same α/τ axes as the
    summary landscape so the figures read together.
    """
    if full.empty:
        return []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt = plot_format.lstrip(".")

    apply_theme()

    df_clean = full.dropna(subset=["Abs_Gap"])
    if df_clean.empty:
        return []
    vmin = float(df_clean["Abs_Gap"].min())
    vmax = float(df_clean["Abs_Gap"].max())

    species_present = [s for s in _SPECIES_ORDER if s in df_clean["Species"].unique()]
    if not species_present:
        species_present = sorted(df_clean["Species"].unique())

    paths: list[Path] = []
    for sp in species_present:
        sp_df = df_clean[df_clean["Species"] == sp]
        if sp_df.empty:
            continue
        genes = sorted(sp_df["Gene"].unique())
        n = len(genes)
        cols = min(4, n)
        rows = math.ceil(n / cols)
        fig, axes = plt.subplots(rows, cols,
                                  figsize=(3.2 * cols + 1.2, 2.9 * rows + 0.8),
                                  squeeze=False)
        for ax in axes.flat[n:]:
            ax.axis("off")

        im = None
        for idx, (ax, gene) in enumerate(zip(axes.flat, genes)):
            sub = sp_df[sp_df["Gene"] == gene]
            pivot = sub.pivot(index="tau", columns="alpha", values="Abs_Gap")
            pivot = pivot.iloc[::-1]
            row, col = divmod(idx, cols)
            im = _draw_heatmap(
                ax, pivot, vmin=vmin, vmax=vmax, title=gene,
                draw_xlabel=(row == rows - 1),
                draw_ylabel=(col == 0),
            )

        fig.suptitle(
            f"{sp} — per-gene CoFold parameter landscape\n"
            "mean |CoFold − DMS| ΔG;  dark = smaller gap = better fit",
            fontsize=TICK_FONTSIZE + 1, y=1.02, fontweight="bold",
        )
        fig.tight_layout(rect=[0, 0, 0.91, 1.0])
        if im is not None:
            cbar_ax = fig.add_axes([0.93, 0.15, 0.018, 0.70])
            cbar = fig.colorbar(im, cax=cbar_ax)
            cbar.set_label("|CoFold − DMS| ΔG  (kcal/mol)", fontsize=TICK_FONTSIZE)
            cbar.ax.tick_params(labelsize=TICK_FONTSIZE)

        out_path = out_dir / f"cofold_parameter_landscape_{sp.lower()}.{fmt}"
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        paths.append(out_path)
    return paths
