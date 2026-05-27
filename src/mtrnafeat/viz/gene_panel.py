"""Per-gene composition / paired-pair / local-foldedness panel.

Generalizes the legacy ND6 special-case panel to every gene in the .db
files. For each (species, gene) we emit a 3-panel figure:
  1. Base counts (A / U / G / C)
  2. Paired-pair composition (G-C / A-U / G-U wobble) as % of paired columns
  3. Local paired-fraction track along the transcript with a dedicated
     architecture strip below it (5'UTR / CDS / 3'UTR) so the region
     labels never overlap the trace title or the data.

The title reports overall sequence GC% so cross-gene comparisons are easy.

Layout note: the third panel is split into two stacked sub-axes via a
nested gridspec — a tall trace axis and a thin architecture strip. This
puts all region labels in their own dedicated band, eliminating the
"5'UTR/3'UTR labels overlap the trace title" issue the previous
shade_regions-above-the-axis approach had.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from mtrnafeat.analysis.statistics import paired_composition, sequence_gc_pct
from mtrnafeat.constants import PALETTE
from mtrnafeat.io.db_parser import DbRecord
from mtrnafeat.viz.style import (
    LABEL_FONTSIZE,
    TITLE_FONTSIZE,
    apply_theme,
    style_axis,
)


def _draw_architecture_strip(ax, l_utr5: int, l_cds: int, transcript_len: int) -> None:
    """Compact architecture strip: colored rectangles only; labels go in figure legend."""
    cds_start = l_utr5 + 1
    cds_end = l_utr5 + l_cds
    utr3_start = cds_end + 1

    regions = []
    if l_utr5 >= 1:
        regions.append((1, l_utr5, PALETTE.get("UTR", "#dddddd")))
    if cds_end >= cds_start:
        regions.append((cds_start, cds_end, PALETTE.get("CDS", "#7faaff")))
    if transcript_len >= utr3_start:
        regions.append((utr3_start, transcript_len, PALETTE.get("UTR", "#dddddd")))

    ax.set_xlim(1, transcript_len)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", left=False, labelleft=False)

    for lo, hi, color in regions:
        width = hi - lo + 1
        ax.add_patch(plt.Rectangle((lo, 0.18), width, 0.64, facecolor=color,
                                    edgecolor="#444444", linewidth=0.6, zorder=1))


def plot_gene(rec: DbRecord, species: str, out_path: Path,
              annot: dict | None = None, dpi: int = 300) -> Path:
    """3-panel per-gene composition + structure figure.

    `annot` (optional) carries `l_utr5`, `l_cds`, `l_tr` from io/annotations.py;
    when present, the third panel gets a dedicated architecture strip below
    the paired-fraction trace.
    """
    apply_theme()
    seq = rec.sequence
    gene = rec.gene
    counts = {n: seq.count(n) for n in "AUGC"}
    gc_p, au_p, gu_p = paired_composition(seq, rec.structure)
    n = len(seq)

    fig = plt.figure(figsize=(16, 5.8), constrained_layout=False)
    outer = fig.add_gridspec(
        1, 3, width_ratios=[1, 1, 2.4], wspace=0.28,
        left=0.05, right=0.98, top=0.90, bottom=0.12,
    )

    ax_counts = fig.add_subplot(outer[0])
    ax_pairs = fig.add_subplot(outer[1])

    has_annot = (
        annot is not None
        and int(annot.get("l_utr5", 0)) + int(annot.get("l_cds", 0)) > 0
    )
    if has_annot:
        right = outer[2].subgridspec(2, 1, height_ratios=[7.5, 1.0], hspace=0.06)
        ax_track = fig.add_subplot(right[0])
        ax_arch = fig.add_subplot(right[1], sharex=ax_track)
    else:
        ax_track = fig.add_subplot(outer[2])
        ax_arch = None

    bars = ax_counts.bar(list(counts.keys()), list(counts.values()),
                          color=["#aec7e8", "#1f77b4", "#ff7f0e", "#d62728"])
    ax_counts.set_title("Base counts", fontsize=13, fontweight="bold", pad=8)
    ax_counts.set_ylabel("Count")
    ax_counts.margins(y=0.08)
    for b, v in zip(bars, counts.values(), strict=True):
        ax_counts.text(b.get_x() + b.get_width() / 2,
                        v + max(1, max(counts.values()) * 0.015),
                        str(v), ha="center", fontsize=10)
    style_axis(ax_counts)

    ax_pairs.bar(["G-C", "A-U", "G-U wobble"], [gc_p, au_p, gu_p],
                  color=["#2ca02c", "#9467bd", "#e377c2"])
    ax_pairs.set_ylabel("% of paired columns")
    ax_pairs.set_title("Paired-pair composition", fontsize=13, fontweight="bold", pad=8)
    ax_pairs.set_ylim(0, 105)
    for i, v in enumerate([gc_p, au_p, gu_p]):
        ax_pairs.text(i, v + 1.5, f"{v:.1f}%", ha="center", fontsize=10)
    style_axis(ax_pairs)

    paired = np.array([1 if ch in "()" else 0 for ch in rec.structure])
    win = 50
    if win >= n:
        win = max(2, n // 4)
    smooth = np.convolve(paired, np.ones(win) / win, mode="same")
    ax_track.plot(range(1, n + 1), smooth, color="black", linewidth=1.6, zorder=3)
    ax_track.fill_between(range(1, n + 1), smooth, 0, alpha=0.3, color="darkblue", zorder=2)
    ax_track.set_ylim(0, 1)
    ax_track.set_xlim(1, n)
    ax_track.set_ylabel("DMS-guided paired fraction", fontsize=LABEL_FONTSIZE - 1)
    ax_track.set_title(f"DMS-guided paired fraction (window = {win} nt)",
                        fontsize=13, fontweight="bold", pad=8)
    style_axis(ax_track)

    if ax_arch is not None:
        ax_track.tick_params(labelbottom=False)
        ax_track.set_xlabel("")
        try:
            _draw_architecture_strip(
                ax_arch,
                l_utr5=int(annot.get("l_utr5", 0)),
                l_cds=int(annot.get("l_cds", 0)),
                transcript_len=int(annot.get("l_tr", n)),
            )
        except Exception:
            pass
        ax_arch.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE - 1)
        utr_patch = mpatches.Patch(facecolor=PALETTE.get("UTR", "#dddddd"),
                                    edgecolor="#444444", linewidth=0.6, label="5'/3' UTR")
        cds_patch = mpatches.Patch(facecolor=PALETTE.get("CDS", "#7faaff"),
                                    edgecolor="#444444", linewidth=0.6, label="CDS")
        ax_track.legend(handles=[utr_patch, cds_patch], loc="upper right",
                         fontsize=9, frameon=True, framealpha=0.85,
                         edgecolor="#cccccc", title="Region", title_fontsize=8)
    else:
        ax_track.set_xlabel("Transcript position (nt)", fontsize=LABEL_FONTSIZE - 1)

    fig.suptitle(f"{species} {gene} — sequence GC = {sequence_gc_pct(seq):.1f}%",
                 fontsize=TITLE_FONTSIZE, y=0.98, fontweight="bold")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


# Backward-compat alias: the previous ND6-specific function.
def plot_nd6(rec: DbRecord, out_path: Path, dpi: int = 300) -> Path:
    return plot_gene(rec, species=getattr(rec, "species", "Sample"),
                      out_path=out_path, annot=None, dpi=dpi)
