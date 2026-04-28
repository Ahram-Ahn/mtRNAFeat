"""Tests for the RNAstructure engine wrapper.

The pure-Python helpers (`validate_rna_alphabet`, `dms_to_minus999`) run
unconditionally. The actual `fold()` test is gated on RNAstructure being
installed (Fold binary on PATH and DATAPATH set); it skips cleanly otherwise.
"""
from __future__ import annotations

import os
import shutil

import numpy as np
import pytest

from mtrnafeat.engines._common import (
    MissingEngineError,
    dms_to_minus999,
    validate_rna_alphabet,
)
from mtrnafeat.engines import rnastructure


def _has_rnastructure() -> bool:
    if "DATAPATH" not in os.environ:
        return False
    if os.environ.get("RNASTRUCTURE_FOLD"):
        return True
    binroot = os.environ.get("RNASTRUCTURE_BIN")
    if binroot and (os.path.exists(os.path.join(binroot, "Fold"))):
        return True
    return shutil.which("Fold") is not None


needs_rnastructure = pytest.mark.skipif(
    not _has_rnastructure(),
    reason="RNAstructure 'Fold' not on PATH or DATAPATH unset",
)


# ---------------------------------------------------------------------------
# Pure-Python helpers
# ---------------------------------------------------------------------------

def test_validate_rna_alphabet_accepts_acgu():
    validate_rna_alphabet("ACGUACGU")


def test_validate_rna_alphabet_rejects_empty():
    with pytest.raises(ValueError, match="Empty"):
        validate_rna_alphabet("")


def test_validate_rna_alphabet_rejects_t():
    with pytest.raises(ValueError, match="T detected"):
        validate_rna_alphabet("ACGT")


def test_validate_rna_alphabet_rejects_lowercase():
    with pytest.raises(ValueError, match="Lowercase"):
        validate_rna_alphabet("acgu")


def test_validate_rna_alphabet_rejects_n():
    with pytest.raises(ValueError, match="Non-RNA"):
        validate_rna_alphabet("ACNG")


def test_dms_to_minus999_masks_g_and_u():
    seq = "ACGU"
    dms = np.array([0.5, 0.6, 0.7, 0.8])
    out = dms_to_minus999(dms, seq)
    assert out[0] == 0.5         # A passes through
    assert out[1] == 0.6         # C passes through
    assert out[2] == -999.0      # G masked
    assert out[3] == -999.0      # U masked


def test_dms_to_minus999_masks_nan():
    seq = "ACAC"
    dms = np.array([0.5, np.nan, 0.7, 0.8])
    out = dms_to_minus999(dms, seq)
    assert out[0] == 0.5
    assert out[1] == -999.0
    assert out[2] == 0.7
    assert out[3] == 0.8


def test_dms_to_minus999_length_check():
    with pytest.raises(ValueError, match="length"):
        dms_to_minus999(np.array([0.1, 0.2]), "ACG")


# ---------------------------------------------------------------------------
# Subprocess wrapper — gated on RNAstructure being installed
# ---------------------------------------------------------------------------

def test_fold_raises_missing_engine_when_unavailable(monkeypatch):
    """If neither RNASTRUCTURE_FOLD nor RNASTRUCTURE_BIN nor PATH find Fold,
    we should get a MissingEngineError with an actionable message."""
    monkeypatch.delenv("RNASTRUCTURE_FOLD", raising=False)
    monkeypatch.delenv("RNASTRUCTURE_BIN", raising=False)
    monkeypatch.setenv("PATH", "")  # nothing on PATH
    with pytest.raises(MissingEngineError, match="Fold"):
        rnastructure.fold("GGGAAACCC")


@needs_rnastructure
def test_fold_returns_balanced_dot_bracket():
    """30 nt of repeated mini-hairpin sequence; assert parens balance."""
    seq = "GGGAAAUCCC" * 3
    res = rnastructure.fold(seq, max_bp_span=30)
    assert isinstance(res.dot_bracket, str)
    assert len(res.dot_bracket) == len(seq)
    assert res.dot_bracket.count("(") == res.dot_bracket.count(")")
    assert all(ch in ".()" for ch in res.dot_bracket)
    # Real folds on this sequence are negative (real hairpins).
    assert res.mfe <= 0.0


@needs_rnastructure
def test_fold_max_bp_span_excludes_long_range_pairs():
    """A short max_bp_span should produce a structure with no long-range pairs."""
    seq = "GGGGGAAAAAAAAAACCCCC" * 3   # 60 nt; long-range pairing possible without -md
    span = 15
    res = rnastructure.fold(seq, max_bp_span=span)
    # parse parens to pair table and check spans
    stack: list[int] = []
    pairs: list[tuple[int, int]] = []
    for i, ch in enumerate(res.dot_bracket):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            pairs.append((j, i))
    assert all(b - a <= span for a, b in pairs), \
        f"found pair exceeding span={span}: {pairs}"
