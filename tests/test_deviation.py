"""Pure-python tests for the structure-deviation analysis layer.

These tests do not invoke ViennaRNA — region calling, merging, and
classification are exercised on synthetic deviation tracks so the
deterministic logic is locked in independently of folding.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mtrnafeat.analysis.deviation import (
    DeviationResult,
    _bin_intervals,
    _empirical_p_one_sided_max,
    _filter_length,
    _merge_runs,
    _runs,
    annotate_regions,
    call_regions,
    classify_region,
    gene_region_matrix,
    gene_summary,
    paired_binary_from_dotbracket,
    per_position_table,
)


# ──────────────────── parser / building blocks ────────────────────


def test_paired_binary_round_trip():
    assert paired_binary_from_dotbracket("((..))").tolist() == [1, 1, 0, 0, 1, 1]


def test_runs_basic():
    mask = np.array([0, 1, 1, 0, 1, 1, 1, 0, 0, 1], dtype=bool)
    assert _runs(mask) == [(1, 3), (4, 7), (9, 10)]
    assert _runs(np.array([], dtype=bool)) == []


def test_merge_runs_respects_gap():
    runs = [(0, 5), (8, 12), (40, 50)]
    # gap 0..4 absorbs (8, 12) into (0, 12); (40, 50) stays separate
    assert _merge_runs(runs, gap=4) == [(0, 12), (40, 50)]
    # gap=2 keeps (0, 5) and (8, 12) separate
    assert _merge_runs(runs, gap=2) == [(0, 5), (8, 12), (40, 50)]


def test_filter_length_drops_short():
    assert _filter_length([(0, 5), (10, 40)], min_len=10) == [(10, 40)]


# ──────────────────── region calling ────────────────────


def test_call_regions_recovers_synthetic_intervals():
    # A 200-nt track with one positive plateau (20:80) and one negative
    # plateau (120:180), both well above the 0.25 threshold.
    n = 200
    dev = np.zeros(n)
    dev[20:80] = 0.45
    dev[120:180] = -0.40
    out = call_regions(dev, threshold=0.25, min_region_length=25, merge_gap=10)
    assert len(out) == 2
    # Sort by start
    out.sort(key=lambda r: r[0])
    (s1, e1, sign1), (s2, e2, sign2) = out
    assert (s1, e1, sign1) == (20, 80, "pos")
    assert (s2, e2, sign2) == (120, 180, "neg")


def test_call_regions_drops_short():
    n = 100
    dev = np.zeros(n)
    dev[40:60] = 0.4   # 20 nt — below min_region_length
    out = call_regions(dev, threshold=0.25, min_region_length=25, merge_gap=10)
    assert out == []


def test_call_regions_merges_adjacent_same_sign():
    n = 200
    dev = np.zeros(n)
    dev[20:50] = 0.4
    dev[55:90] = 0.4   # 5-nt gap; merge_gap=10 should fuse them
    out = call_regions(dev, threshold=0.25, min_region_length=25, merge_gap=10)
    # One merged region
    assert [(s, e) for s, e, _ in out] == [(20, 90)]


# ──────────────────── classification ────────────────────


def test_classify_region_all_branches():
    high = 0.50
    low = 0.30
    dev = 0.25
    assert classify_region(0.7, 0.1, +0.6, high=high, low=low,
                           dev_thresh=dev)[0] == "model_high_dms_low"
    assert classify_region(0.1, 0.7, -0.6, high=high, low=low,
                           dev_thresh=dev)[0] == "model_low_dms_high"
    # Strong deviation but model & dms don't both pass thresholds → mixed
    assert classify_region(0.55, 0.20, +0.35, high=high, low=low,
                           dev_thresh=dev)[0] == "model_high_dms_low"
    assert classify_region(0.50, 0.40, +0.30, high=high, low=low,
                           dev_thresh=dev)[0] == "mixed_deviation"
    # Weak deviation, both high → concordant_paired
    assert classify_region(0.7, 0.7, +0.0, high=high, low=low,
                           dev_thresh=dev)[0] == "concordant_paired"
    # Weak deviation, both low → concordant_open
    assert classify_region(0.1, 0.1, +0.0, high=high, low=low,
                           dev_thresh=dev)[0] == "concordant_open"
    # Weak deviation in the middle → ambiguous
    assert classify_region(0.40, 0.40, +0.0, high=high, low=low,
                           dev_thresh=dev)[0] == "ambiguous"


# ──────────────────── annotation + summary ────────────────────


def _make_result(p_model, p_dms, *, sequence=None, structure=None,
                 species="Test", gene="X", smooth=5):
    p_model = np.asarray(p_model, dtype=float)
    p_dms = np.asarray(p_dms, dtype=float)
    n = p_model.size
    if sequence is None:
        sequence = "A" * n
    if structure is None:
        structure = "." * n
    pms = pd.Series(p_model).rolling(smooth, center=True, min_periods=1).mean().to_numpy()
    pds = pd.Series(p_dms).rolling(smooth, center=True, min_periods=1).mean().to_numpy()
    return DeviationResult(
        species=species, gene=gene, sequence=sequence, dms_structure=structure,
        p_model_raw=p_model, p_model_smooth=pms,
        p_dms_raw=p_dms, p_dms_smooth=pds,
        deviation_raw=p_model - p_dms,
        deviation_smooth=pms - pds,
        rolling_window=smooth,
        rnaplfold_window=80, rnaplfold_max_bp_span=50, rnaplfold_cutoff=0.001,
    )


def test_per_position_table_distance_and_region():
    n = 60
    p = np.full(n, 0.5)
    res = _make_result(p, p)
    annot = {"l_tr": n, "l_utr5": 10, "l_cds": 30, "l_utr3": 20}
    df = per_position_table(res, annot)
    # First base is upstream of CDS by 10 nt
    assert df.loc[0, "Transcript_Region"] == "5UTR"
    assert df.loc[0, "Distance_To_Start"] == -10  # cds_start_1 = 11
    # Position 11 (1-based) is the first CDS base, distance 0
    cds_first = df[df["Distance_To_Start"] == 0].iloc[0]
    assert cds_first["Transcript_Region"] == "CDS"
    assert cds_first["CDS_Phase"] == 0


def test_annotate_regions_classifies_and_overlaps_tis():
    """A synthetic transcript with the start codon at position 50 (0-based)
    and an open patch around it should produce a model_high_dms_low region
    that overlaps the TIS window."""
    n = 200
    p_model = np.full(n, 0.85)
    p_dms = np.full(n, 0.85)
    # Open the TIS: model stays high, DMS drops
    p_dms[40:65] = 0.05
    res = _make_result(p_model, p_dms, smooth=5)
    annot = {"l_tr": n, "l_utr5": 50, "l_cds": 120, "l_utr3": 30}
    pp = per_position_table(res, annot)
    transcript_region = pp["Transcript_Region"].tolist()
    regions = call_regions(res.deviation_smooth, threshold=0.25,
                           min_region_length=10, merge_gap=10)
    assert regions, "expected at least one region"

    class _Cfg:
        structure_deviation_high_threshold = 0.50
        structure_deviation_low_threshold = 0.30
        structure_deviation_threshold = 0.25
        structure_deviation_tis_upstream = 30
        structure_deviation_tis_downstream = 60
        structure_deviation_stop_upstream = 60
        structure_deviation_stop_downstream = 30

    df = annotate_regions(res, regions, annot, transcript_region, cfg=_Cfg())
    assert (df["Region_Class"] == "model_high_dms_low").any()
    tis_overlapping = df[df["Overlaps_TIS_Window"]]
    assert not tis_overlapping.empty, "expected the called region to overlap the TIS window"


def test_bin_intervals_partitions_short_transcripts_safely():
    """Short transcripts should still produce non-negative-width bins
    or empty bins — never (lo, hi) with hi < lo."""
    n = 100

    class _Cfg:
        structure_deviation_tis_upstream = 30
        structure_deviation_tis_downstream = 60
        structure_deviation_stop_upstream = 60
        structure_deviation_stop_downstream = 30
        structure_deviation_early_cds_nt = 300

    annot = {"l_tr": n, "l_utr5": 0, "l_cds": 90, "l_utr3": 10}
    bins = _bin_intervals(annot, n, cfg=_Cfg())
    for name, (lo, hi) in bins.items():
        assert lo >= 0
        assert hi <= n
        assert hi >= 0  # may equal lo (empty bin), but never negative width
        assert lo <= hi or lo == hi


def test_gene_summary_counts_regions_by_class():
    n = 300
    p_model = np.full(n, 0.85)
    p_dms = np.full(n, 0.85)
    p_dms[100:160] = 0.05  # one strong open patch
    res = _make_result(p_model, p_dms, smooth=5)
    annot = {"l_tr": n, "l_utr5": 50, "l_cds": 200, "l_utr3": 50}

    class _Cfg:
        structure_deviation_high_threshold = 0.50
        structure_deviation_low_threshold = 0.30
        structure_deviation_threshold = 0.25
        structure_deviation_tis_upstream = 30
        structure_deviation_tis_downstream = 60
        structure_deviation_stop_upstream = 60
        structure_deviation_stop_downstream = 30
        structure_deviation_early_cds_nt = 300

    pp = per_position_table(res, annot)
    regions = call_regions(res.deviation_smooth, threshold=0.25,
                           min_region_length=10, merge_gap=10)
    reg_df = annotate_regions(res, regions, annot,
                              pp["Transcript_Region"].tolist(), cfg=_Cfg())
    summary = gene_summary(res, reg_df, annot, cfg=_Cfg())
    assert summary["N_Regions_Total"] >= 1
    assert summary["N_Model_High_DMS_Low"] == int(
        (reg_df["Region_Class"] == "model_high_dms_low").sum()
    )
    assert summary["TIS_Mean_Model_Paired"] >= 0.0


def test_empirical_p_one_sided_max_matches_definition():
    # Null distribution {0.1, 0.2, 0.3, 0.4, 0.5}. Observed = 0.35 → 2 of 5
    # null max-stats are >= observed, so p = (2 + 1) / (5 + 1) = 0.5.
    null = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    p = _empirical_p_one_sided_max(0.35, null)
    assert p == pytest.approx(3 / 6)


def test_empirical_p_one_sided_max_observed_above_all():
    null = np.array([0.1, 0.2, 0.3])
    # Observed strictly larger than every null value → p = 1/(N+1).
    p = _empirical_p_one_sided_max(0.99, null)
    assert p == pytest.approx(1 / 4)


def test_empirical_p_one_sided_max_handles_nan_inputs():
    # NaN observed → NaN p; all-NaN null → NaN p; empty null → NaN p.
    assert np.isnan(_empirical_p_one_sided_max(float("nan"), np.array([0.1])))
    assert np.isnan(_empirical_p_one_sided_max(0.5, np.array([np.nan, np.nan])))
    assert np.isnan(_empirical_p_one_sided_max(0.5, np.array([], dtype=float)))


def test_annotate_regions_with_null_populates_p_and_label():
    """A region whose observed |dev| sits below most of a null distribution
    should be flagged ``max_stat_nonsignificant``; the inverse should be
    flagged ``max_stat_significant``."""
    n = 100
    p_model = np.full(n, 0.85)
    p_dms = np.full(n, 0.85)
    # One strong open patch around positions 30-60 → big |dev|
    p_dms[30:60] = 0.05
    res = _make_result(p_model, p_dms, smooth=5)
    annot = {"l_tr": n, "l_utr5": 20, "l_cds": 60, "l_utr3": 20}

    class _Cfg:
        structure_deviation_high_threshold = 0.50
        structure_deviation_low_threshold = 0.30
        structure_deviation_threshold = 0.25
        structure_deviation_tis_upstream = 30
        structure_deviation_tis_downstream = 60
        structure_deviation_stop_upstream = 60
        structure_deviation_stop_downstream = 30
        structure_deviation_null_model = "dinuc"

    pp = per_position_table(res, annot)
    regions = call_regions(res.deviation_smooth, threshold=0.25,
                           min_region_length=10, merge_gap=10)
    assert regions, "expected at least one region in the open patch"

    # Null where every value sits below the observed |dev| and there are
    # enough draws that 1/(N+1) drops below the 0.05 significance cutoff.
    null_small = np.full(100, 0.05)
    df_sig = annotate_regions(
        res, regions, annot, pp["Transcript_Region"].tolist(),
        cfg=_Cfg(), null_max_stats=null_small,
    )
    assert (df_sig["Statistical_Support_Label"] == "max_stat_significant").all()
    assert (df_sig["Empirical_P"] < 0.05).all()
    assert (df_sig["Null_Model"] == "dinuc").all()
    assert (df_sig["N_Null"] == 100).all()

    # Huge null → observed |dev| is well within the null range → ns
    null_large = np.full(50, 5.0)
    df_ns = annotate_regions(
        res, regions, annot, pp["Transcript_Region"].tolist(),
        cfg=_Cfg(), null_max_stats=null_large,
    )
    assert (df_ns["Statistical_Support_Label"] == "max_stat_nonsignificant").all()


def test_annotate_regions_without_null_keeps_effect_size_only_label():
    n = 100
    p_model = np.full(n, 0.85)
    p_dms = np.full(n, 0.85)
    p_dms[30:60] = 0.05
    res = _make_result(p_model, p_dms, smooth=5)
    annot = {"l_tr": n, "l_utr5": 20, "l_cds": 60, "l_utr3": 20}

    class _Cfg:
        structure_deviation_high_threshold = 0.50
        structure_deviation_low_threshold = 0.30
        structure_deviation_threshold = 0.25
        structure_deviation_tis_upstream = 30
        structure_deviation_tis_downstream = 60
        structure_deviation_stop_upstream = 60
        structure_deviation_stop_downstream = 30

    pp = per_position_table(res, annot)
    regions = call_regions(res.deviation_smooth, threshold=0.25,
                           min_region_length=10, merge_gap=10)
    df = annotate_regions(
        res, regions, annot, pp["Transcript_Region"].tolist(), cfg=_Cfg(),
    )
    assert (df["Statistical_Support_Label"] == "effect_size_only").all()
    assert df["Empirical_P"].isna().all()
    assert (df["N_Null"] == 0).all()


def test_gene_region_matrix_rows_present():
    n = 200
    p = np.full(n, 0.5)
    res = _make_result(p, p)
    annot = {"l_tr": n, "l_utr5": 30, "l_cds": 150, "l_utr3": 20}

    class _Cfg:
        structure_deviation_tis_upstream = 30
        structure_deviation_tis_downstream = 60
        structure_deviation_stop_upstream = 60
        structure_deviation_stop_downstream = 30
        structure_deviation_early_cds_nt = 300

    df = gene_region_matrix(res, annot, cfg=_Cfg())
    # At minimum TIS / mid_CDS / stop_proximal exist; 5_end and 3_end may
    # be empty when the CDS dominates.
    assert "TIS" in df["Region_Bin"].tolist()
    # Every reported bin has at least one position
    assert (df["N_Positions"] > 0).all()
