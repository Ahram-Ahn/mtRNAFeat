"""Single source of plotting style and helpers. Publication-ready defaults."""
from __future__ import annotations

import math
from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from mtrnafeat.constants import PALETTE


def plot_path(out_dir, basename: str, fmt: str = "png") -> Path:
    """Compose `<out_dir>/<basename>.<fmt>` so commands don't hard-code an extension.

    matplotlib infers the backend (PNG/SVG/PDF) from the file extension, so
    this is the single switch point for changing the export format across the
    whole package via cfg.plot_format.
    """
    ext = (fmt or "png").lstrip(".").lower()
    return Path(out_dir) / f"{basename}.{ext}"


SPINE_WIDTH = 1.2
TICK_WIDTH = 1.1
TICK_LENGTH = 4.5
TITLE_FONTSIZE = 17
LABEL_FONTSIZE = 14
TICK_FONTSIZE = 12
LEGEND_FONTSIZE = 12
LINEWIDTH = 2.4
AXIS_COLOR = "#222222"


def apply_theme() -> None:
    """Publication context. Call once at the top of every plot function."""
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "axes.titlesize": TITLE_FONTSIZE,
        "axes.labelsize": LABEL_FONTSIZE,
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "xtick.labelsize": TICK_FONTSIZE,
        "ytick.labelsize": TICK_FONTSIZE,
        "legend.fontsize": LEGEND_FONTSIZE,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def style_axis(ax) -> None:
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE_WIDTH)
        spine.set_color(AXIS_COLOR)
    ax.tick_params(axis="both", which="major", width=TICK_WIDTH,
                   length=TICK_LENGTH, labelsize=TICK_FONTSIZE)


def shade_regions(data_axes, l_utr5: int, l_cds: int, transcript_len: int,
                    label_axis=None, fontsize: int = 11) -> None:
    """Tint UTR/CDS/3'UTR backgrounds across one or more shared-x data panels.

    Replaces the old `add_region_track` strip (which made short 5'UTRs
    illegible) with: subtle alpha-tinted axvspans on every passed `data_axes`
    + dashed vertical boundary lines + region labels placed above the top
    panel. For very short regions (<5% of transcript length) the label is
    drawn outside the band with a thin leader line pointing into it.

    Pass `label_axis` to control where labels appear (defaults to the first
    axis in `data_axes`). Pass a list of axes for `data_axes` if you have
    multiple stacked panels and want them all tinted consistently.
    """
    if not isinstance(data_axes, (list, tuple)):
        data_axes = [data_axes]
    if label_axis is None:
        label_axis = data_axes[0]

    cds_start = l_utr5 + 1
    cds_end = l_utr5 + l_cds
    utr3_start = cds_end + 1

    regions = []
    if l_utr5 >= 1:
        regions.append(("5'UTR", 1, l_utr5, PALETTE.get("UTR", "#dddddd")))
    if cds_end >= cds_start:
        regions.append(("CDS", cds_start, cds_end, PALETTE.get("CDS", "#7faaff")))
    if transcript_len >= utr3_start:
        regions.append(("3'UTR", utr3_start, transcript_len, PALETTE.get("UTR", "#dddddd")))

    for name, lo, hi, color in regions:
        for ax in data_axes:
            ax.axvspan(lo, hi + 0.5, color=color, alpha=0.18 if name != "CDS" else 0.10, zorder=0)

    for boundary in (l_utr5 + 0.5, cds_end + 0.5):
        if 0 < boundary < transcript_len:
            for ax in data_axes:
                ax.axvline(boundary, color="black", linestyle="--", linewidth=0.7, alpha=0.55,
                            zorder=1)

    # Place labels above the top panel.
    label_axis.set_xlim(1, transcript_len)
    short_threshold = max(40, int(0.05 * transcript_len))
    for name, lo, hi, _ in regions:
        center = (lo + hi) / 2.0
        width = hi - lo + 1
        # Use axes-fraction y so label sits just above the top of the panel.
        if width >= short_threshold:
            label_axis.annotate(
                name, xy=(center, 1.0), xycoords=("data", "axes fraction"),
                xytext=(0, 6), textcoords="offset points",
                ha="center", va="bottom", fontsize=fontsize, fontweight="bold",
                color="#333333",
            )
        else:
            # Short region — leader line from above into the band.
            leader_y_off = 22 + (8 if name.startswith("5") else 0)
            label_axis.annotate(
                name, xy=(center, 1.0), xycoords=("data", "axes fraction"),
                xytext=(0, leader_y_off), textcoords="offset points",
                ha="center", va="bottom", fontsize=fontsize - 1, fontweight="bold",
                color="#333333",
                arrowprops=dict(arrowstyle="-", color="#444444", lw=0.6,
                                  shrinkA=0, shrinkB=2),
            )


def add_region_track(ax, l_utr5: int, l_cds: int, transcript_len: int,
                      x_min: float | None = None, x_max: float | None = None) -> None:
    if x_min is None:
        x_min = 1
    if x_max is None:
        x_max = transcript_len
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, 1)
    ax.axis("off")
    cds_start = l_utr5 + 1
    cds_end = l_utr5 + l_cds
    utr3_start = cds_end + 1
    if l_utr5 >= 1:
        ax.add_patch(Rectangle((1, 0.18), l_utr5, 0.64, color=PALETTE["UTR"], ec=None))
        ax.text((1 + l_utr5) / 2, 0.5, "5'UTR", ha="center", va="center",
                fontsize=11, color="#444444")
    if cds_end >= cds_start:
        ax.add_patch(Rectangle((cds_start, 0.18), cds_end - cds_start + 1, 0.64,
                                color=PALETTE["CDS"], ec=None))
        ax.text((cds_start + cds_end) / 2, 0.5, "CDS", ha="center", va="center",
                fontsize=12, fontweight="bold", color="white")
    if transcript_len >= utr3_start:
        ax.add_patch(Rectangle((utr3_start, 0.18), transcript_len - utr3_start + 1, 0.64,
                                color=PALETTE["UTR"], ec=None))
        ax.text((utr3_start + transcript_len) / 2, 0.5, "3'UTR", ha="center", va="center",
                fontsize=11, color="#444444")


def repel_labels(ax, xs: Iterable[float], ys: Iterable[float], labels: Iterable[str],
                  *, color: str = "black", fontsize: int = 10,
                  k: int = 24, max_iter: int = 400, halo: bool = True,
                  use_adjusttext: bool = True) -> None:
    """Force-based label de-collision with white halo + leader lines.

    If `adjustText` is importable and `use_adjusttext` is True we delegate to
    it (it does collision detection against scatter points + lines + iso-
    contours, much better than this pure-matplotlib fallback). Otherwise we
    use the in-house spring solver.

    Halo: applied via `path_effects` (white stroke under the text) so labels
    stay legible when drawn over dense KDE clouds.

    `k` is the spring constant in display-pixel units; lower → looser packing.
    """
    import matplotlib.patheffects as path_effects

    xs = list(xs); ys = list(ys); labels = list(labels)
    if not xs:
        return

    halo_effects = ([path_effects.Stroke(linewidth=3.2, foreground="white"),
                      path_effects.Normal()] if halo else None)

    if use_adjusttext:
        try:
            from adjustText import adjust_text  # type: ignore
            text_objs = []
            for x, y, lbl in zip(xs, ys, labels):
                t = ax.text(x, y, lbl, fontsize=fontsize, color=color, fontweight="bold",
                             ha="center", va="center", zorder=12)
                if halo_effects is not None:
                    t.set_path_effects(halo_effects)
                text_objs.append(t)
            adjust_text(text_objs, ax=ax,
                          arrowprops=dict(arrowstyle="-", color=color, lw=0.6, alpha=0.6),
                          expand_points=(1.4, 1.4), expand_text=(1.2, 1.2))
            return
        except ImportError:
            pass

    fig = ax.figure
    fig.canvas.draw()
    inv = ax.transData.inverted()

    pts = ax.transData.transform(list(zip(xs, ys)))
    label_pts = [tuple(p) for p in pts]
    radii = []
    text_objs = []
    for (x, y), lbl in zip(label_pts, labels):
        t = ax.text(0, 0, lbl, fontsize=fontsize, color=color, fontweight="bold",
                     ha="center", va="center", zorder=12)
        if halo_effects is not None:
            t.set_path_effects(halo_effects)
        bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
        radii.append(0.5 * math.hypot(bb.width, bb.height) + 4.0)
        text_objs.append(t)

    moved = list(label_pts)
    for _ in range(max_iter):
        any_overlap = False
        for i in range(len(moved)):
            xi, yi = moved[i]
            ri = radii[i]
            fx = fy = 0.0
            for j in range(len(moved)):
                if i == j:
                    continue
                dx = xi - moved[j][0]
                dy = yi - moved[j][1]
                dist = math.hypot(dx, dy) or 1e-3
                min_d = ri + radii[j]
                if dist < min_d:
                    any_overlap = True
                    push = (min_d - dist) / dist
                    fx += dx * push * 0.5
                    fy += dy * push * 0.5
            ax_dx = xi - pts[i][0]
            ax_dy = yi - pts[i][1]
            d_anchor = math.hypot(ax_dx, ax_dy) or 1e-3
            anchor_pull = max(0.0, d_anchor - 18.0) / 140.0  # gentler pull
            fx -= ax_dx * anchor_pull
            fy -= ax_dy * anchor_pull
            moved[i] = (xi + fx, yi + fy)
        if not any_overlap:
            break

    for (lx, ly), (px, py), t in zip(moved, pts, text_objs):
        data_xy = inv.transform((lx, ly))
        t.set_position((data_xy[0], data_xy[1]))
        if math.hypot(lx - px, ly - py) > 6.0:
            anchor_data = inv.transform((px, py))
            ax.plot([anchor_data[0], data_xy[0]], [anchor_data[1], data_xy[1]],
                    color=color, linewidth=0.6, alpha=0.6, zorder=11)


def species_legend_handles(species: list[str]) -> list[Line2D]:
    return [Line2D([0], [0], marker="o", color="w",
                    markerfacecolor=PALETTE.get(sp, "black"),
                    markersize=11, markeredgecolor="black", label=sp)
            for sp in species]
