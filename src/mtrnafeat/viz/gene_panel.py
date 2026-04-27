"""Per-gene composition / paired-pair / local-foldedness panel.

Generalizes the legacy ND6 special-case panel to every gene in the .db
files. For each (species, gene) we emit a 3-panel figure:
  1. Base counts (A / U / G / C)
  2. Paired-pair composition (G-C / A-U / G-U wobble) as % of paired columns
  3. Local paired-fraction track along the transcript with a region track
     (5'UTR / CDS / 3'UTR) overlay where annotations are available.

Title reports overall sequence GC% so cross-gene comparisons are easy.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mtrnafeat.analysis.statistics import paired_composition, sequence_gc_pct
from mtrnafeat.io.db_parser import DbRecord
from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    TITLE_FONTSIZE,
    apply_theme,
    shade_regions,
    style_axis,
)


def plot_gene(rec: DbRecord, species: str, out_path: Path,
              annot: dict | None = None, dpi: int = 300) -> Path:
    """3-panel per-gene composition + structure figure.

    `annot` (optional) carries `l_utr5`, `l_cds`, `l_tr` from io/annotations.py;
    when present, the third panel gets a region-track overlay.
    """
    apply_theme()
    seq = rec.sequence
    gene = rec.gene
    counts = {n: seq.count(n) for n in "AUGC"}
    gc_p, au_p, gu_p = paired_composition(seq, rec.structure)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5),
                              gridspec_kw={"width_ratios": [1, 1, 2.4]})

    bars = axes[0].bar(list(counts.keys()), list(counts.values()),
                        color=["#aec7e8", "#1f77b4", "#ff7f0e", "#d62728"])
    axes[0].set_title("Base counts", fontsize=12, fontweight="bold")
    axes[0].set_ylabel("Count")
    for b, v in zip(bars, counts.values()):
        axes[0].text(b.get_x() + b.get_width() / 2, v + max(1, max(counts.values()) * 0.01),
                      str(v), ha="center", fontsize=10)
    style_axis(axes[0])

    axes[1].bar(["G-C", "A-U", "G-U wobble"], [gc_p, au_p, gu_p],
                 color=["#2ca02c", "#9467bd", "#e377c2"])
    axes[1].set_ylabel("% of paired columns")
    axes[1].set_title("Paired-pair composition", fontsize=12, fontweight="bold")
    axes[1].set_ylim(0, 100)
    for i, v in enumerate([gc_p, au_p, gu_p]):
        axes[1].text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=10)
    style_axis(axes[1])

    n = len(seq)
    paired = np.array([1 if ch in "()" else 0 for ch in rec.structure])
    win = max(20, n // 30)
    if win >= n:
        win = max(2, n // 4)
    smooth = np.convolve(paired, np.ones(win) / win, mode="same")
    axes[2].plot(range(1, n + 1), smooth, color="black", linewidth=1.6, zorder=3)
    axes[2].fill_between(range(1, n + 1), smooth, 0, alpha=0.3, color="darkblue", zorder=2)
    axes[2].set_ylim(0, 1)
    axes[2].set_xlim(1, n)
    axes[2].set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE)
    axes[2].set_ylabel("Local paired fraction")
    axes[2].set_title(f"Local paired fraction (window={win} nt)",
                       fontsize=12, fontweight="bold")
    style_axis(axes[2])

    if annot is not None:
        try:
            shade_regions([axes[2]], l_utr5=int(annot.get("l_utr5", 0)),
                          l_cds=int(annot.get("l_cds", 0)),
                          transcript_len=int(annot.get("l_tr", n)),
                          label_axis=axes[2])
        except Exception:
            pass

    fig.suptitle(f"{species} {gene} — sequence GC = {sequence_gc_pct(seq):.1f}%",
                 fontsize=TITLE_FONTSIZE, y=1.03)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)


# Backward-compat alias: the previous ND6-specific function.
def plot_nd6(rec: DbRecord, out_path: Path, dpi: int = 300) -> Path:
    return plot_gene(rec, species=getattr(rec, "species", "Sample"),
                      out_path=out_path, annot=None, dpi=dpi)
