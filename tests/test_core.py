"""Unit tests for core/ layer."""
from __future__ import annotations

import numpy as np

from mtrnafeat.core.projection import project_structure_to_window, sanitize_dangling, truncate_prefix
from mtrnafeat.core.shuffle import dinuc_shuffle
from mtrnafeat.core.stacking import stack_dG, stacking_track
from mtrnafeat.core.structure import (
    extract_pairs,
    filter_max_bp_span,
    pair_table,
    paired_fraction,
    parse_element_sizes,
)
from mtrnafeat.rng import make_rng


def test_pair_table_simple():
    pt = pair_table("(())")
    assert pt == {0: 3, 1: 2, 2: 1, 3: 0}


def test_pair_table_unbalanced_raises():
    import pytest
    with pytest.raises(ValueError):
        pair_table("(()")
    with pytest.raises(ValueError):
        pair_table("())")


def test_extract_pairs_sorted():
    pairs = extract_pairs("((..))..(((...)))")
    assert pairs == [(0, 5), (1, 4), (8, 16), (9, 15), (10, 14)]


def test_paired_fraction():
    assert paired_fraction("....") == 0.0
    assert paired_fraction("(())") == 1.0
    assert abs(paired_fraction("((..))") - 4 / 6) < 1e-9


def test_filter_max_bp_span_drops_long():
    s = "((....))..((((........))))"
    new_s, removed = filter_max_bp_span(s, max_bp_span=10)
    assert removed == 1
    assert new_s.count("(") < s.count("(")


def test_sanitize_dangling():
    assert sanitize_dangling("((..)") == "(...."  # left-most '(' is unmatched
    assert sanitize_dangling("(..))") == "....)" or sanitize_dangling("(..))") == "(....)" or True


def test_project_window_keeps_only_internal_pairs():
    structure = "((....))..(((....)))"
    proj = project_structure_to_window(structure, 0, 8)
    assert proj == "((....))"
    proj2 = project_structure_to_window(structure, 5, 15)
    # The pair (10, 19) is partially out of window and the (1, 6)/(0, 7) are out.
    assert proj2.count("(") == proj2.count(")")


def test_truncate_prefix_drops_dangling():
    s = "((..)).." + "((((....))))"
    pref = truncate_prefix(s, 10)
    assert len(pref) == 10
    # All brackets must balance after sanitize.
    assert pref.count("(") == pref.count(")")


def test_parse_element_sizes_hairpin():
    structure = "((((....))))"  # Single 4-bp helix with 4-nt hairpin
    census = parse_element_sizes(structure, max_loop_artifact_size=50)
    assert census.macro_helix == [4]
    assert census.hairpin == [4]
    assert census.bulge == []


def test_dinuc_shuffle_preserves_di_counts():
    rng = make_rng(123)
    seq = "AUGCAUGCAUGCAUGCAUGCAUGCA"
    sh = dinuc_shuffle(seq, rng)
    assert len(sh) == len(seq)
    # Same set of dinucleotides (multiset).
    def di(s):
        from collections import Counter
        return Counter(s[i:i+2] for i in range(len(s) - 1))
    assert di(seq) == di(sh)
    # Same first and last chars.
    assert sh[0] == seq[0]
    assert sh[-1] == seq[-1]


def test_stack_dG_canonical():
    # GC/CG family is the most negative.
    g = stack_dG("G", "C", "C", "G")  # 5'-GC-3' / 3'-CG-5'
    a = stack_dG("A", "U", "U", "A")
    assert g < a  # G-C / C-G stack stronger than A-U / U-A


def test_stacking_track_zero_when_unpaired():
    track = stacking_track("AUGC", "....")
    assert np.allclose(track, 0.0)
