"""Tests that exercise the analysis layer without ViennaRNA."""
from __future__ import annotations

import pytest

from mtrnafeat.analysis.statistics import paired_composition, sequence_gc_pct, transcript_stats
from mtrnafeat.io.db_parser import parse_db


def test_sequence_gc_pct():
    assert sequence_gc_pct("AAGCAA") == pytest.approx(2 / 6 * 100)
    assert sequence_gc_pct("") == 0.0


def test_paired_composition_pair_typed():
    seq = "GCAUGCAU"
    struct = "((((..))" + "))"  # not a valid balanced string
    # Use an obviously balanced one:
    seq2 = "GCAUUUUGC"
    struct2 = "(((...)))"
    gc, au, gu = paired_composition(seq2, struct2)
    # Pairs from outside-in: (G, C) -> GC, (C, G) -> GC, (A, U) -> AU
    assert gc == pytest.approx(200.0 / 3)  # 2/3
    assert au == pytest.approx(100.0 / 3)
    assert gu == 0.0


def test_transcript_stats_columns(mini_human_db):
    rec = parse_db(mini_human_db)[0]
    row = transcript_stats(rec, condition="X")
    expected = {"Condition", "Gene", "Length", "MFE", "Normalized_MFE_per_nt",
                "Foldedness_Pct", "Sequence_GC_Pct", "Paired_GC_Pct",
                "Paired_AU_Pct", "Paired_GU_Pct"}
    assert expected.issubset(row.keys())
    assert row["Length"] == 538
    assert row["MFE"] == -26.3
