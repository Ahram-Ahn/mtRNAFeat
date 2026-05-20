"""CoFold parameter-sweep figures.

Three outputs (all species on one figure each):

1. ``cofold_gap_closure`` — the core narrative figure.  For each species
   (Human | Yeast, side-by-side panels), plots the *fraction of the
   Vienna→DMS ΔG gap that is closed* as alpha increases, with one line
   per τ and ±1 SD shading across genes.  A horizontal reference at
   frac=1.0 marks exact agreement with DMS.  Human lines stay below 1
   even at α=1 (needs stronger penalty); Yeast lines cross 1 at α≈0.75
   (moderate penalty suffices).

2. ``cofold_parameter_landscape`` — alpha × tau heatmap, color = mean
   |CoFold − DMS| across genes, one panel per species side-by-side.
   The cell with the smallest mean gap is marked with a star so the
   optimum is immediately obvious.

3. ``cofold_per_window_rmse`` — per-gene window-level RMSE curves
   (α on x, RMSE on y, one line per τ), replacing the Pearson-r curves
   which barely move with the penalty.  Lower RMSE = the local ΔG
   profile shape + magnitude both track DMS more closely.
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
    """Add 'frac_closed' column: (CoFold_MFE − Vienna_α0) / (DMS_Eval_dG − Vienna_α0).

    0 = no improvement over plain Vienna, 1 = exactly matches DMS ΔG.
    Values > 1 mean the penalty over-corrects (too folded relative to DMS).
    """
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
# Figure 1: gap-closure convergence (the main narrative)
# ---------------------------------------------------------------------------

def _draw_closure_panel(ax, species_df: pd.DataFrame, species: str,
                         taus: list, cmap: list) -> None:
    """Draw gap-closure curves for one species onto `ax`."""
    agg = (
        species_df.groupby(["alpha", "tau"])["frac_closed"]
        .agg(["mean", "std"])
        .reset_index()
        .rename(columns={"mean": "mu", "std": "sd"})
    )
    alphas_sorted = sorted(agg["alpha"].unique())

    # Overshoot shading (anything above frac=1 is over-penalised)
    ax.axhspan(1.0, 1.6, alpha=0.07, color="#CC0000", zorder=0, lw=0)

    for color, tau in zip(cmap, taus):
        sub = agg[agg["tau"] == tau].sort_values("alpha")
        mu = sub["mu"].to_numpy()
        sd = sub["sd"].fillna(0).to_numpy()
        xs = sub["alpha"].to_numpy()
        ax.plot(xs, mu, "-o", color=color, lw=2.0, ms=6,
                label=f"τ = {int(tau)} nt")
        ax.fill_between(xs, mu - sd, mu + sd, alpha=0.13, color=color)

    # DMS target line
    ax.axhline(1.0, color="#333333", ls="--", lw=1.3, zorder=5)
    ax.text(
        max(alphas_sorted) * 0.98, 1.035,
        "CoFold = DMS ΔG",
        ha="right", va="bottom", fontsize=TICK_FONTSIZE,
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
    """Side-by-side gap-closure convergence figure (Human | Yeast).

    This is the primary narrative figure: shows at what alpha the
    CoFold MFE matches the DMS-evaluated ΔG, revealing that human
    requires a stronger long-range penalty than yeast.
    """
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

    # Shared τ legend on the rightmost panel
    handles, labels = axes[0][-1].get_legend_handles_labels()
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

def gap_heatmap_panels(full: pd.DataFrame, out_dir: Path, plot_format: str,
                        dpi: int = 300) -> list[Path]:
    """Side-by-side α × τ heatmaps of mean |CoFold − DMS| (Human | Yeast).

    Dark cells = small gap (good fit).  The optimum cell is marked with ★.
    Shows at a glance that human requires (high α, low τ) while yeast
    optimum sits at a lower α.
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

    # Shared colour scale across species for fair comparison
    pivot_all = (
        full.dropna(subset=["Abs_Gap"])
        .groupby(["Species", "alpha", "tau"])["Abs_Gap"]
        .mean()
        .reset_index()
    )
    vmin = pivot_all["Abs_Gap"].min()
    vmax = pivot_all["Abs_Gap"].max()

    for ax, sp, letter in zip(axes[0], species_present, "AB"):
        sub = pivot_all[pivot_all["Species"] == sp]
        pivot = sub.pivot(index="tau", columns="alpha", values="Abs_Gap")
        # Reverse tau axis so small tau (= early penalty) is at the top.
        pivot = pivot.iloc[::-1]

        im = ax.imshow(
            pivot.values,
            aspect="auto",
            cmap="YlOrRd_r",   # reversed: dark = low gap = good
            vmin=vmin, vmax=vmax,
            origin="upper",
        )

        alphas = [f"{a:.2g}" for a in pivot.columns.tolist()]
        taus = [str(int(t)) for t in pivot.index.tolist()]
        ax.set_xticks(range(len(alphas)))
        ax.set_xticklabels(alphas, fontsize=TICK_FONTSIZE)
        ax.set_yticks(range(len(taus)))
        ax.set_yticklabels(taus, fontsize=TICK_FONTSIZE)
        ax.set_xlabel("α — penalty strength (kcal/mol)", fontsize=LABEL_FONTSIZE)
        ax.set_ylabel("τ — decay length (nt)", fontsize=LABEL_FONTSIZE)
        ax.set_title(sp, fontsize=TITLE_FONTSIZE, fontweight="bold", pad=6)
        panel_label(ax, letter)

        # Annotate each cell with the mean gap value
        for r in range(pivot.shape[0]):
            for c in range(pivot.shape[1]):
                val = pivot.values[r, c]
                if np.isfinite(val):
                    ax.text(c, r, f"{val:.1f}", ha="center", va="center",
                            fontsize=6.5, color="white" if val < vmax * 0.55 else "#222222")

        # Mark optimum cell
        flat_idx = np.nanargmin(pivot.values)
        best_r, best_c = np.unravel_index(flat_idx, pivot.shape)
        ax.text(best_c, best_r, "★", ha="center", va="center",
                fontsize=14, color="#FFD700", fontweight="bold",
                path_effects=_star_outline())

    fig.suptitle(
        "CoFold parameter landscape — mean gap to DMS-evaluated ΔG\n"
        "★ = optimal (α, τ) per species;  dark = smaller gap = better fit",
        fontsize=TICK_FONTSIZE, y=1.02,
    )
    fig.tight_layout(rect=[0, 0, 0.88, 1.0])

    # Colorbar in the reserved right margin
    cbar_ax = fig.add_axes([0.90, 0.15, 0.025, 0.70])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Mean |CoFold − DMS| ΔG  (kcal/mol)", fontsize=TICK_FONTSIZE)
    cbar.ax.tick_params(labelsize=TICK_FONTSIZE)

    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return [out_path]


def _star_outline():
    import matplotlib.patheffects as pe
    return [pe.Stroke(linewidth=2.0, foreground="#333333"), pe.Normal()]


# Backwards-compatible alias (cofold command used to call gap_strip_panels).
gap_strip_panels = gap_heatmap_panels


# ---------------------------------------------------------------------------
# Figure 3: per-window RMSE curves (replaces Pearson-r grid)
# ---------------------------------------------------------------------------

def _rmse_grid_for_species(species_win: pd.DataFrame, species: str,
                             out_path: Path, dpi: int) -> Path:
    apply_theme()
    if species_win.empty:
        return out_path
    genes = sorted(species_win["Gene"].unique())
    n = len(genes)
    cols = min(4, n)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(5.8 * cols, 3.2 * rows), squeeze=False)
    for ax in axes.flat[n:]:
        ax.axis("off")

    taus = sorted(species_win["tau"].unique())
    cmap = sns.color_palette("cividis", len(taus))

    for ax, gene in zip(axes.flat, genes):
        sub = species_win[species_win["Gene"] == gene]
        for color, tau in zip(cmap, taus):
            line = sub[sub["tau"] == tau].sort_values("alpha")
            ax.plot(line["alpha"], line["Per_Window_RMSE"],
                    "-o", color=color, lw=1.6, ms=4.5, label=f"τ={int(tau)}")
        ax.set_xlabel("α", fontsize=LABEL_FONTSIZE - 2)
        ax.set_ylabel("Window RMSE (kcal/mol)", fontsize=LABEL_FONTSIZE - 2)
        ax.set_title(gene, fontsize=TITLE_FONTSIZE - 3, fontweight="bold")
        ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.4)
        ax.set_axisbelow(True)
        style_axis(ax)
        legend_outside(ax, position="right", fontsize=7, frameon=False,
                       title="τ (nt)", title_fontsize=7)

    fig.suptitle(
        f"{species} — per-window RMSE: CoFold ΔG vs. DMS-projected ΔG across α sweep\n"
        "Lower RMSE = local energy profile magnitude tracks DMS more closely",
        fontsize=TICK_FONTSIZE + 1, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def per_window_corr_curves(win: pd.DataFrame, out_dir: Path, plot_format: str,
                             dpi: int = 300) -> list[Path]:
    """Per-gene window-level RMSE curves (α on x), one figure per species."""
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
        out_path = out_dir / f"cofold_per_window_rmse_{sp.lower()}.{fmt}"
        _rmse_grid_for_species(win[win["Species"] == sp], sp, out_path, dpi)
        paths.append(out_path)
    return paths
