"""TIS zoom visualization: −50 / +50 nt around the start codon.

Two-panel layout: Human on the left, Yeast on the right.

Genes whose 5'UTR is shorter than `upstream_nt` get a hatched bar so the
reader knows the window's upstream context is truncated. The legend
shows the actual upstream length used per gene via the x-tick labels:
`GENE (5'UTR=N)`.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from mtrnafeat.viz.style import apply_theme, style_axis, LABEL_FONTSIZE, TITLE_FONTSIZE


_SPECIES_ORDER = ["Human", "Yeast"]


def _draw_one(ax, sub: pd.DataFrame, species: str) -> None:
    if sub.empty:
        ax.text(0.5, 0.5, f"No {species} data.", ha="center", va="center",
                transform=ax.transAxes)
        ax.axis("off")
        return
    sub = sub.copy().reset_index(drop=True)
    ind = range(len(sub))
    width = 0.42
    full_ctx = sub["Has_Full_5UTR_Context"].astype(bool).values \
        if "Has_Full_5UTR_Context" in sub.columns else [True] * len(sub)
    dms_bars = ax.bar([i - width/2 for i in ind], sub["DMS_TIS_Energy"], width=width,
                       color="#1F77B4", label="DMS-projected ΔG")
    mfe_bars = ax.bar([i + width/2 for i in ind], sub["Vienna_TIS_MFE"], width=width,
                       color="#D62728", label="Vienna MFE")
    for bars in (dms_bars, mfe_bars):
        for bar, full in zip(bars, full_ctx):
            if not full:
                bar.set_hatch("///")
                bar.set_edgecolor("white")
    # Tick labels include actual upstream context length.
    if "L_5UTR_in_window" in sub.columns:
        labels = [f"{g}\n(5'UTR={int(u)})" for g, u in zip(sub["Gene"], sub["L_5UTR_in_window"])]
    else:
        labels = list(sub["Gene"])
    ax.set_xticks(list(ind))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_title(f"{species} TIS (−50 / +50 nt around start)",
                 fontsize=TITLE_FONTSIZE - 2, pad=10)
    ax.set_ylabel("Energy (kcal/mol)", fontsize=LABEL_FONTSIZE)
    ax.legend(fontsize=10, loc="best")
    style_axis(ax)


def tis_zoom_panel(df_tis: pd.DataFrame, out_path: Path, dpi: int = 300) -> Path:
    apply_theme()
    if df_tis.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No TIS data.", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)
    df = df_tis.dropna(subset=["DMS_TIS_Energy", "Vienna_TIS_MFE"]).copy()
    species_present = [s for s in _SPECIES_ORDER if s in df["Species"].unique()]
    if not species_present:
        species_present = sorted(df["Species"].unique())
    n = len(species_present)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 5.5), sharey=False)
    if n == 1:
        axes = [axes]
    for ax, sp in zip(axes, species_present):
        _draw_one(ax, df[df["Species"] == sp].sort_values("Gene"), sp)
    fig.suptitle("TIS zoom — DMS-projected ΔG vs Vienna MFE\n"
                  "(hatched bars = 5'UTR shorter than 50 nt → upstream context truncated)",
                  fontsize=TITLE_FONTSIZE - 1, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
