"""Small statistical helpers shared across analysis stages.

Kept dependency-free (pure numpy) so we don't pull in statsmodels just
for a handful of multiple-testing corrections.
"""
from __future__ import annotations

import numpy as np


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg step-up FDR adjustment.

    Returns q-values aligned with the input. NaN p-values pass through
    as NaN and are excluded from the rank denominator (so a partial
    table doesn't penalize the tested rows). Result is clipped to [0, 1].
    """
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    pf = p[finite]
    n = pf.size
    order = np.argsort(pf, kind="mergesort")
    ranked = pf[order]
    # Step-up: q_(i) = min over k>=i of (n/k) * p_(k)
    raw = ranked * n / np.arange(1, n + 1)
    q_sorted = np.minimum.accumulate(raw[::-1])[::-1]
    q = np.empty(n, dtype=float)
    q[order] = np.clip(q_sorted, 0.0, 1.0)
    out[finite] = q
    return out
