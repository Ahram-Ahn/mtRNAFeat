"""Sample-label helpers for plots that accept arbitrary db_files labels."""
from __future__ import annotations

from collections.abc import Iterable

import seaborn as sns

from mtrnafeat.constants import PALETTE

_FALLBACK_PALETTE = (
    "#2166AC", "#D6604D", "#1B9E77", "#7570B3",
    "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
)


def ordered_samples(values: Iterable[str]) -> list[str]:
    """Unique sample labels, keeping Human/Yeast first when present."""
    seen: list[str] = []
    for value in values:
        label = str(value)
        if label not in seen:
            seen.append(label)
    preferred = [s for s in ("Human", "Yeast") if s in seen]
    return preferred + [s for s in seen if s not in preferred and s != "n/a"]


def sample_palette(samples: Iterable[str]) -> dict[str, str]:
    labels = ordered_samples(samples)
    colors = sns.color_palette(_FALLBACK_PALETTE, max(len(labels), 1)).as_hex()
    return {
        label: PALETTE.get(label, colors[i % len(colors)])
        for i, label in enumerate(labels)
    }


def sample_color(sample: str, fallback_index: int = 0) -> str:
    if sample in PALETTE:
        return PALETTE[sample]
    return _FALLBACK_PALETTE[fallback_index % len(_FALLBACK_PALETTE)]
