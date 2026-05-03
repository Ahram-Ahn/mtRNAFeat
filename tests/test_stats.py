"""Tests for mtrnafeat.core.stats helpers."""
from __future__ import annotations

import numpy as np

from mtrnafeat.core.stats import bh_fdr


def test_bh_fdr_matches_textbook_example():
    # Classic BH worked example: p = [0.01, 0.02, 0.03, 0.04, 0.05]
    # at n=5, the q-values are p_i * 5 / rank, then step-up min from the
    # right. This collapses to [0.05, 0.05, 0.05, 0.05, 0.05].
    p = np.array([0.01, 0.02, 0.03, 0.04, 0.05])
    q = bh_fdr(p)
    np.testing.assert_allclose(q, [0.05, 0.05, 0.05, 0.05, 0.05], atol=1e-12)


def test_bh_fdr_preserves_input_order():
    p = np.array([0.05, 0.01, 0.04, 0.02, 0.03])
    q = bh_fdr(p)
    # The smallest p (0.01 at index 1) maps to the smallest q.
    assert np.argmin(q) == 1


def test_bh_fdr_handles_nan_passthrough():
    p = np.array([0.01, np.nan, 0.02, 0.05, np.nan])
    q = bh_fdr(p)
    assert np.isnan(q[1]) and np.isnan(q[4])
    # The three finite p-values should be rank-corrected as a 3-element family.
    finite_q = q[[0, 2, 3]]
    assert np.all(np.isfinite(finite_q))
    assert np.all((finite_q >= 0) & (finite_q <= 1))


def test_bh_fdr_clipped_to_unit_interval():
    # All p-values are 1.0 — every q should be exactly 1.0 (clipped).
    p = np.ones(5)
    q = bh_fdr(p)
    np.testing.assert_allclose(q, np.ones(5))


def test_bh_fdr_all_nan_returns_all_nan():
    p = np.full(4, np.nan)
    q = bh_fdr(p)
    assert np.all(np.isnan(q))


def test_bh_fdr_monotone_in_sorted_p():
    # When p-values are sorted ascending, q-values must be non-decreasing
    # (BH is a step-up procedure on the sorted ranks).
    rng = np.random.default_rng(0)
    p = np.sort(rng.uniform(size=20))
    q = bh_fdr(p)
    assert np.all(np.diff(q) >= -1e-12)
