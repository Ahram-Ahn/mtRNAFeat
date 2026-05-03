"""DrTransformer kinetic-trajectory plot: occupancy vs transcript length."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

from mtrnafeat.viz.style import LABEL_FONTSIZE, TITLE_FONTSIZE, apply_theme, style_axis


def plot_kinetic_trajectory(traj_df, gene: str, species: str, out_path: Path,
                              top_n: int = 8, dpi: int = 300) -> Path:
    apply_theme()
    sub = traj_df[(traj_df["Gene"] == gene) & (traj_df["Species"] == species)].copy()
    if sub.empty:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, f"No kinetic trajectory for {species} {gene}.",
                ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=dpi)
        plt.close(fig)
        return Path(out_path)
    final_L = sub["Transcript_Length"].max()
    top_ids = (sub[sub["Transcript_Length"] == final_L]
               .sort_values("Occupancy", ascending=False)["Struct_ID"].head(top_n).tolist())
    fig, axes = plt.subplots(2, 1, figsize=(13, 8.5), sharex=True)
    palette = sns.color_palette("colorblind", n_colors=max(top_n, 4))
    for c, sid in enumerate(top_ids):
        s = sub[sub["Struct_ID"] == sid].sort_values("Transcript_Length")
        axes[0].plot(s["Transcript_Length"], s["Occupancy"], color=palette[c], lw=1.8, label=sid)
        axes[1].plot(s["Transcript_Length"], s["Energy"],    color=palette[c], lw=1.8, label=sid)
    axes[0].set_ylabel("Occupancy")
    axes[0].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=10)
    axes[0].set_title(f"{species} {gene}: kinetic trajectory (top {top_n} structures)",
                      fontsize=TITLE_FONTSIZE, pad=12)
    axes[1].set_ylabel("Energy (kcal/mol)")
    axes[1].set_xlabel("Transcript length (nt)", fontsize=LABEL_FONTSIZE)
    for a in axes:
        style_axis(a)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return Path(out_path)
