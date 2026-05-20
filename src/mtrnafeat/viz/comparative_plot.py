"""Yeast↔human COX1 comparative plots: substitution bar chart + directional flux."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from mtrnafeat.viz.style import LABEL_FONTSIZE, TICK_FONTSIZE, TITLE_FONTSIZE, apply_theme


def plot_directional_flux(flux_df, out_path: Path, dpi: int = 300) -> Path:
    """12 (from→to) × 3 (codon position) heatmap of bias score.

    Rows sorted by |mean bias| across positions (most biased on top).
    BH-significant cells outlined in black (q<0.05).
    """
    apply_theme()
    if flux_df.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No directional flux data.", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)
    df = flux_df.copy()
    df["Direction"] = df["From"] + "→" + df["To"]
    pivot = df.pivot_table(index="Direction", columns="Position", values="Bias_Score", aggfunc="mean")
    qpivot = df.pivot_table(index="Direction", columns="Position", values="BH_q", aggfunc="mean")

    # Sort rows by absolute mean bias across all positions (most biased on top)
    row_order = pivot.abs().mean(axis=1).sort_values(ascending=False).index
    pivot = pivot.loc[row_order]
    qpivot = qpivot.reindex(row_order)

    fig, ax = plt.subplots(figsize=(7.5, 8))
    sns.heatmap(pivot, ax=ax, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                annot_kws={"size": 10}, linewidths=0.4,
                cbar_kws={"label": "Bias score (obs − exp)"})

    # Outline BH-significant cells
    for i, direction in enumerate(pivot.index):
        for j, pos in enumerate(pivot.columns):
            try:
                q = qpivot.loc[direction, pos]
            except KeyError:
                q = 1.0
            if q < 0.05:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                                            edgecolor="black", lw=2.0))

    ax.set_title(
        "Yeast→Human substitution flux (BH-significant outlined, q<0.05)\n"
        "Positive = enriched;  negative = depleted vs null expectation",
        fontsize=TITLE_FONTSIZE - 2, pad=10,
    )
    ax.set_xlabel("Codon position  (1 = first,  2 = second,  3 = wobble)",
                  fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Substitution direction", fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


def plot_substitution_summary(sub_summary, out_path: Path, dpi: int = 300) -> Path:
    """12-state substitution spectrum (Yeast → Human COX1) as 2×3 heatmaps.

    Rows = Identical AA / Divergent AA. Columns = codon position 1 / 2 / 3.
    Each tile is a 4×4 (Yeast base × Human base) heatmap. Diagonal cells
    are unchanged nucleotides and labelled "Same" — excluded from the
    12-state count. Off-diagonal cells show count + percent-within-tile.
    """
    apply_theme()
    if sub_summary.empty:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "No substitutions observed.", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)

    bases = ["A", "C", "G", "T"]
    df = sub_summary.copy()
    # Normalize U to T
    df["Yeast_Base"] = df["Yeast_Base"].replace({"U": "T"})
    df["Human_Base"] = df["Human_Base"].replace({"U": "T"})

    df["AA_Cat"] = df["Same_AA"].map({True: "Identical AA", False: "Divergent AA"})
    aa_order = ["Identical AA", "Divergent AA"]
    positions = [1, 2, 3]

    # Build the per-tile percentage normalization across off-diagonal cells.
    tile_totals: dict[tuple[str, int], float] = {}
    for cat in aa_order:
        for pos in positions:
            sub = df[(df["AA_Cat"] == cat) & (df["Position"] == pos)]
            off = sub[sub["Yeast_Base"] != sub["Human_Base"]]
            tile_totals[(cat, pos)] = float(off["Count"].sum())

    # Vmax for the colorbar uses the maximum off-diagonal percentage.
    all_pcts: list[float] = []
    for cat in aa_order:
        for pos in positions:
            sub = df[(df["AA_Cat"] == cat) & (df["Position"] == pos)]
            tot = tile_totals[(cat, pos)] or 1.0
            for _, r in sub.iterrows():
                if r["Yeast_Base"] != r["Human_Base"]:
                    all_pcts.append(100.0 * float(r["Count"]) / tot)
    vmax = max(all_pcts) if all_pcts else 1.0
    if vmax <= 0:
        vmax = 1.0

    import matplotlib as mpl
    cmap = plt.get_cmap("viridis")
    norm = mpl.colors.Normalize(vmin=0.0, vmax=vmax)

    fig, axes = plt.subplots(
        2, 3, figsize=(15, 9.5),
        gridspec_kw={"hspace": 0.32, "wspace": 0.18},
    )

    for ri, cat in enumerate(aa_order):
        for ci, pos in enumerate(positions):
            ax = axes[ri][ci]
            sub = df[(df["AA_Cat"] == cat) & (df["Position"] == pos)]
            tot = tile_totals[(cat, pos)] or 1.0

            counts = np.zeros((4, 4), dtype=int)
            for _, r in sub.iterrows():
                yi = bases.index(r["Yeast_Base"]) if r["Yeast_Base"] in bases else None
                hi = bases.index(r["Human_Base"]) if r["Human_Base"] in bases else None
                if yi is None or hi is None:
                    continue
                counts[yi, hi] += int(r["Count"])

            pct = np.zeros((4, 4), dtype=float)
            for yi in range(4):
                for hi in range(4):
                    if yi == hi:
                        continue
                    pct[yi, hi] = 100.0 * counts[yi, hi] / tot

            # Draw cells manually to allow "Same" labels on the diagonal.
            for yi in range(4):
                for hi in range(4):
                    if yi == hi:
                        rect_color = "#EEEEEE"
                        ax.add_patch(plt.Rectangle((hi - 0.5, yi - 0.5), 1, 1,
                                                    facecolor=rect_color,
                                                    edgecolor="white", lw=1.0))
                        ax.text(hi, yi, "Same", ha="center", va="center",
                                fontsize=9, color="#999999", style="italic")
                        continue
                    c = cmap(norm(pct[yi, hi]))
                    ax.add_patch(plt.Rectangle((hi - 0.5, yi - 0.5), 1, 1,
                                                facecolor=c, edgecolor="white", lw=1.0))
                    # Text color depends on luminance of fill
                    lum = 0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]
                    text_color = "white" if lum < 0.55 else "#1a1a1a"
                    ax.text(hi, yi, f"{counts[yi, hi]}\n({pct[yi, hi]:.1f}%)",
                            ha="center", va="center", fontsize=9,
                            fontweight="bold", color=text_color)

            ax.set_xticks(range(4))
            ax.set_yticks(range(4))
            ax.set_xticklabels(bases, fontsize=TICK_FONTSIZE + 1, fontweight="bold")
            ax.set_yticklabels(bases, fontsize=TICK_FONTSIZE + 1, fontweight="bold")
            ax.set_xlim(-0.5, 3.5)
            ax.set_ylim(3.5, -0.5)
            ax.set_aspect("equal")
            ax.set_title(f"{cat} | Position {pos}",
                          fontsize=TITLE_FONTSIZE - 1, fontweight="bold", pad=6)
            if ri == 1:
                ax.set_xlabel("Human nucleotide", fontsize=LABEL_FONTSIZE)
            if ci == 0:
                ax.set_ylabel("Yeast nucleotide", fontsize=LABEL_FONTSIZE)
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.tick_params(axis="both", length=0)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.93, 0.18, 0.018, 0.64])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Substitution fraction within codon position (%)",
                   fontsize=LABEL_FONTSIZE)
    cbar.ax.tick_params(labelsize=TICK_FONTSIZE)

    fig.suptitle(
        "12-state substitution spectrum in mt-COX1 (Yeast → Human)",
        fontsize=TITLE_FONTSIZE + 2, fontweight="bold", y=0.99,
    )
    fig.text(0.5, 0.05,
             "Diagonal cells mark unchanged nucleotides and are excluded "
             "from the 12-state count.",
             ha="center", fontsize=TICK_FONTSIZE, style="italic",
             color="#555555")
    fig.tight_layout(rect=[0, 0.06, 0.91, 0.96])
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)
