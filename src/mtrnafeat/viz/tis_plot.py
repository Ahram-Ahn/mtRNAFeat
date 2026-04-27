"""5'-end zoom visualization: first 50 nt, DMS vs Vienna MFE per gene.

Two-panel layout: Human on the left, Yeast on the right.

Genes without an annotated 5'UTR (Has_5UTR = False) get a hatched bar so
the reader knows the 50 nt window is pure CDS rather than the textbook
'5'UTR + initiation context'.
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
    has_utr = sub["Has_5UTR"].astype(bool).values if "Has_5UTR" in sub.columns else [True] * len(sub)
    dms_bars = ax.bar([i - width/2 for i in ind], sub["DMS_5end_Energy"], width=width,
                       color="#1F77B4", label="DMS-projected ΔG")
    mfe_bars = ax.bar([i + width/2 for i in ind], sub["Vienna_5end_MFE"], width=width,
                       color="#D62728", label="Vienna MFE")
    for bars in (dms_bars, mfe_bars):
        for bar, has in zip(bars, has_utr):
            if not has:
                bar.set_hatch("///")
                bar.set_edgecolor("white")
    ax.set_xticks(list(ind))
    ax.set_xticklabels(sub["Gene"], rotation=35, ha="right", fontsize=10)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_title(f"{species} — first 50 nt of mRNA",
                 fontsize=TITLE_FONTSIZE - 2, pad=10)
    ax.set_ylabel("Energy (kcal/mol)", fontsize=LABEL_FONTSIZE)
    ax.legend(fontsize=10, loc="best")
    style_axis(ax)


def tis_zoom_panel(df_tis: pd.DataFrame, out_path: Path, dpi: int = 300) -> Path:
    apply_theme()
    if df_tis.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No 5'-end data.", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)
    df = df_tis.dropna(subset=["DMS_5end_Energy", "Vienna_5end_MFE"]).copy()
    species_present = [s for s in _SPECIES_ORDER if s in df["Species"].unique()]
    if not species_present:
        species_present = sorted(df["Species"].unique())
    n = len(species_present)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 5.5), sharey=False)
    if n == 1:
        axes = [axes]
    for ax, sp in zip(axes, species_present):
        _draw_one(ax, df[df["Species"] == sp].sort_values("Gene"), sp)
    fig.suptitle("First 50 nt of mRNA — DMS-projected ΔG vs Vienna MFE\n"
                  "(hatched bars = 5'UTR = 0 → window is pure CDS)",
                  fontsize=TITLE_FONTSIZE - 1, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
