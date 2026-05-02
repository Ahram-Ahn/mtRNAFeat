"""Pure-python tests for the local-probability analysis helpers and the
significance/cotrans provenance columns. None of these tests need
ViennaRNA — they exercise the post-fold reasoning by constructing
``LocalProbResult`` instances directly.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mtrnafeat.analysis import local_probability as lp
from mtrnafeat.analysis.cotrans import add_z_columns, gene_signals  # noqa: F401
from mtrnafeat.analysis.local_probability import (
    LocalProbResult,
    _classify_concordance,
    paired_binary_from_dotbracket,
    per_position_table,
    per_window_agreement_table,
    tis_summary_row,
)


def test_paired_binary_from_dotbracket_basic():
    arr = paired_binary_from_dotbracket("(((...)))")
    assert arr.tolist() == [1, 1, 1, 0, 0, 0, 1, 1, 1]
    assert arr.dtype == np.int8


def test_paired_binary_from_dotbracket_all_unpaired():
    assert paired_binary_from_dotbracket(".....").tolist() == [0, 0, 0, 0, 0]


def test_paired_binary_from_dotbracket_handles_empty():
    arr = paired_binary_from_dotbracket("")
    assert arr.size == 0


def _make_result(p_paired, structure=None, *, species="Test", gene="X"):
    p = np.asarray(p_paired, dtype=float)
    n = p.size
    return LocalProbResult(
        species=species,
        gene=gene,
        sequence="A" * n,
        p_paired=p,
        window=80,
        max_bp_span=50,
        cutoff=0.001,
        dms_structure=structure if structure and len(structure) == n else None,
        dms_paired_binary=(paired_binary_from_dotbracket(structure)
                           if structure and len(structure) == n else None),
    )


def test_per_position_table_includes_dms_columns_when_available():
    result = _make_result([0.1, 0.2, 0.9, 0.95, 0.9, 0.1, 0.2, 0.1, 0.05],
                          structure="(((...)))")
    df = per_position_table(result, annot=None, smooth=3,
                            tis_upstream=10, tis_downstream=10)
    for col in ("P_Paired", "P_Paired_Smoothed",
                "DMS_Paired_Binary", "DMS_Paired_Smoothed",
                "Delta_Ppaired_minus_DMS_Smoothed", "In_TIS_Window"):
        assert col in df.columns
    # DMS binary should match the dot-bracket exactly
    assert df["DMS_Paired_Binary"].astype(int).tolist() == [1, 1, 1, 0, 0, 0, 1, 1, 1]


def test_per_position_table_dms_columns_nan_when_missing():
    # Sequence length 5, no structure attached
    result = _make_result([0.1, 0.5, 0.5, 0.5, 0.1])
    df = per_position_table(result, annot=None)
    assert df["DMS_Paired_Binary"].isna().all()
    assert df["DMS_Paired_Smoothed"].isna().all()
    assert df["Delta_Ppaired_minus_DMS_Smoothed"].isna().all()


def test_per_window_agreement_table_columns_and_values():
    # 12-nt synthetic with a perfect 3-bp stem at each end and an open
    # middle. Make the RNAplfold track agree with the DMS binary.
    structure = "(((......)))"
    p_paired = np.array(
        [0.95, 0.95, 0.95, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.95, 0.95, 0.95]
    )
    result = _make_result(p_paired, structure=structure)
    df = per_window_agreement_table(result, annot=None, win=4, step_nt=4)
    expected_cols = {
        "Species", "Gene", "Window_Start_1based", "Window_End_1based",
        "Window_Center_1based", "Window_Size", "Region_Majority",
        "RNAplfold_Mean_Ppaired", "RNAplfold_Median_Ppaired",
        "RNAplfold_Frac_Low_Ppaired_0p25",
        "DMS_Paired_Fraction", "DMS_Mean_Pair_Span",
        "DMS_Long_Pair_Fraction_gt_50",
        "Agreement_Signed_Delta", "Agreement_Abs_Delta",
    }
    assert expected_cols.issubset(df.columns)
    # First window: positions 1..4 → 3 paired + 1 open. RNAplfold mean
    # should be near (0.95*3 + 0.05) / 4 = 0.725; DMS paired fraction = 0.75.
    first = df.iloc[0]
    assert first["RNAplfold_Mean_Ppaired"] == pytest.approx(0.725, abs=1e-9)
    assert first["DMS_Paired_Fraction"] == pytest.approx(0.75, abs=1e-9)
    # Middle window 5..8 → all open
    mid = df.iloc[1]
    assert mid["DMS_Paired_Fraction"] == pytest.approx(0.0, abs=1e-9)
    assert mid["Agreement_Abs_Delta"] == pytest.approx(0.05, abs=1e-9)


def test_tis_summary_truncates_for_short_5utr():
    """A gene with l_utr5 < tis_upstream must report Has_Full_Upstream=False."""
    n = 200
    p = np.linspace(0.1, 0.9, n)
    structure = "." * n
    result = _make_result(p, structure=structure)
    annot = {"l_tr": n, "l_utr5": 5, "l_cds": 150, "l_utr3": 45}
    rng = np.random.default_rng(7)
    row = tis_summary_row(
        result, annot,
        upstream=50, downstream=50,
        n_circ_shifts=50, rng=rng,
    )
    assert row["TIS_Upstream_Nt_Requested"] == 50
    assert row["TIS_Upstream_Nt_Available"] == 5
    assert row["Has_Full_Upstream_Context"] is False
    # Window starts clamped at the 5' end (1-based)
    assert row["TIS_Window_Start_1based"] == 1
    # Empirical p-values are in (0, 1] (the +1/(n+1) Laplace correction
    # ensures they're never exactly 0)
    p_emp = row["Empirical_P_TIS_Low_Ppaired_vs_CDS_Windows"]
    assert 0.0 < p_emp <= 1.0


def test_tis_summary_full_context_flag_when_5utr_long_enough():
    n = 400
    p = np.full(n, 0.5)
    structure = "." * n
    result = _make_result(p, structure=structure)
    annot = {"l_tr": n, "l_utr5": 100, "l_cds": 250, "l_utr3": 50}
    rng = np.random.default_rng(11)
    row = tis_summary_row(
        result, annot,
        upstream=50, downstream=50,
        n_circ_shifts=20, rng=rng,
    )
    assert row["TIS_Upstream_Nt_Available"] == 50
    assert row["Has_Full_Upstream_Context"] is True


def test_concordance_classes():
    assert _classify_concordance(0.10, 0.10) == "concordant_open"
    assert _classify_concordance(0.80, 0.80) == "concordant_paired"
    assert _classify_concordance(0.10, 0.80) == "RNAplfold_open_only"
    assert _classify_concordance(0.80, 0.10) == "DMS_open_only"
    # 0.30 vs 0.30 is intermediate (neither < 0.25 nor > 0.5)
    assert _classify_concordance(0.30, 0.30) == "discordant"
    assert _classify_concordance(np.nan, 0.10) == "unknown"


def test_cotrans_attach_provenance_tags_every_row():
    from mtrnafeat.analysis.cotrans import attach_provenance

    df = pd.DataFrame({
        "Species": ["Yeast", "Yeast", "Human"],
        "Gene": ["COX1", "COX1", "ND1"],
        "Z_Delta_MFE_per_nt_Smooth": [0.0, -2.5, 1.7],
    })
    tagged = attach_provenance(df)
    assert (tagged["Z_score_type"]
            == "within_gene_window_standardized_delta").all()
    assert (tagged["Is_statistical_pvalue"] == False).all()  # noqa: E712
    # Original frame untouched
    assert "Z_score_type" not in df.columns


def test_cotrans_attach_provenance_handles_empty():
    from mtrnafeat.analysis.cotrans import attach_provenance
    out = attach_provenance(pd.DataFrame())
    assert out.empty


def test_circular_shift_p_low_constant_track_is_one():
    """A constant track means every shift gives the same window mean,
    so the empirical P(shifted ≤ observed) should be 1.0."""
    track = np.full(100, 0.5)
    rng = np.random.default_rng(3)
    p = lp._circular_shift_p_low(track, lo=10, hi=20, n_shifts=50, rng=rng)
    assert p == pytest.approx(1.0)
