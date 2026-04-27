"""Thin tqdm wrapper. Used by every long-running analysis module so that
each step prints a per-task progress bar with ETA and rate.

If tqdm is not installed for some reason, we fall back to a no-op iterator
so the package still functions (just without progress display).
"""
from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import TypeVar

T = TypeVar("T")

try:
    from tqdm.auto import tqdm as _tqdm
    HAS_TQDM = True
except ImportError:  # pragma: no cover
    HAS_TQDM = False
    _tqdm = None


def progress(it: Iterable[T], desc: str | None = None, total: int | None = None,
             leave: bool = True, unit: str = "it", disable: bool = False) -> Iterator[T]:
    """Wrap an iterable with a tqdm progress bar (or pass through if unavailable)."""
    if disable or not HAS_TQDM:
        return iter(it)
    return _tqdm(it, desc=desc, total=total, leave=leave, unit=unit, dynamic_ncols=True)


def step(label: str) -> None:
    """One-shot progress message (printed at start of a stage)."""
    print(f"[mtrnafeat] {label}", flush=True)
