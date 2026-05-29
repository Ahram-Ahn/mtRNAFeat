"""Tests for the redesigned cofold visualisation (no ViennaRNA required).

Uses synthetic DataFrames that mirror the structure of cofold_grid.csv and
cofold_per_window_corr.csv to verify figure output and the core narrative:
human requires a stronger long-range penalty than yeast to match DMS ΔG.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mtrnafeat.viz import cofold_plot

# ---------------------------------------------------------------------------
# Shared synthetic fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def mock_full() -> pd.DataFrame:
    """Synthetic cofold_grid DataFrame: 2 species × 2 genes × 3α × 2τ = 24 rows."""
    rows = []
    for species, dms_frac in [("Human", 0.70), ("Yeast", 0.90)]:
        vienna_mfe = -100.0
        dms_eval = vienna_mfe * dms_frac
        for gene in ["GeneA", "GeneB"]:
            for alpha in [0.0, 0.5, 1.0]:
                for tau in [160.0, 640.0]:
                    penalty = alpha * (1.0 - np.exp(-500.0 / tau)) * 30.0
                    mfe = vienna_mfe + penalty
                    rows.append(dict(
                        Species=species, Gene=gene,
                        alpha=alpha, tau=tau,
                        CoFold_MFE=mfe,
                        DMS_Eval_dG=dms_eval,
                        Gap_CoFold_minus_DMS=mfe - dms_eval,
                        Abs_Gap=abs(mfe - dms_eval),
                    ))
    return pd.DataFrame(rows)


@pytest.fixture(scope="module")
def mock_win() -> pd.DataFrame:
    """Synthetic cofold_per_window_corr DataFrame."""
    rows = []
    for species in ["Human", "Yeast"]:
        for gene in ["GeneA", "GeneB"]:
            for alpha in [0.0, 0.5, 1.0]:
                for tau in [160.0, 640.0]:
                    rows.append(dict(
                        Species=species, Gene=gene,
                        alpha=alpha, tau=tau,
                        Pearson_r_CoFold_vs_DMS=0.70 + alpha * 0.05,
                        Per_Window_RMSE=max(10.0 - alpha * 3.0 - (1 if tau == 160 else 0), 1.0),
                        Per_Window_Mean_Gap=-2.0 + alpha,
                        N_Windows_Valid=10,
                    ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure output tests
# ---------------------------------------------------------------------------

def test_gap_closure_panels_creates_one_file(mock_full, tmp_path):
    paths = cofold_plot.gap_closure_panels(mock_full, tmp_path, "png", dpi=72)
    assert len(paths) == 1
    p = paths[0]
    assert p.exists(), f"File not created: {p}"
    assert p.stat().st_size > 5_000, f"File suspiciously small: {p.stat().st_size} B"
    assert p.name == "cofold_gap_closure.png"


def test_gap_heatmap_panels_creates_one_file(mock_full, tmp_path):
    paths = cofold_plot.gap_heatmap_panels(mock_full, tmp_path, "png", dpi=72)
    assert len(paths) == 1
    p = paths[0]
    assert p.exists()
    assert p.stat().st_size > 5_000
    assert p.name == "cofold_parameter_landscape.png"


def test_per_gene_landscape_creates_per_species_files(mock_full, tmp_path):
    paths = cofold_plot.per_gene_landscape(mock_full, tmp_path, "png", dpi=72)
    assert len(paths) == 2
    names = {p.name for p in paths}
    assert "cofold_parameter_landscape_human.png" in names
    assert "cofold_parameter_landscape_yeast.png" in names
    for p in paths:
        assert p.stat().st_size > 1_000


def test_per_gene_landscape_handles_incomplete_grids(tmp_path):
    """A 5-gene species creates 8 axes; unused axes must not trip strict zip."""
    rows = []
    for gene in ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"]:
        for alpha in [0.0, 0.5]:
            for tau in [160.0, 640.0]:
                rows.append(dict(
                    Species="Yeast", Gene=gene,
                    alpha=alpha, tau=tau,
                    CoFold_MFE=-100 + alpha,
                    DMS_Eval_dG=-90,
                    Gap_CoFold_minus_DMS=-10 + alpha,
                    Abs_Gap=abs(-10 + alpha),
                ))
    paths = cofold_plot.per_gene_landscape(pd.DataFrame(rows), tmp_path, "png", dpi=72)
    assert len(paths) == 1
    assert paths[0].exists()
    assert paths[0].stat().st_size > 1_000


# ---------------------------------------------------------------------------
# Narrative sanity: human needs more penalty than yeast
# ---------------------------------------------------------------------------

def test_narrative_human_needs_stronger_penalty(mock_full):
    """At α=1, yeast should have closed more of the Vienna→DMS gap than human."""
    df = cofold_plot._add_frac_closed(mock_full)
    human = df[(df["Species"] == "Human") & (df["alpha"] == 1.0)]["frac_closed"].mean()
    yeast = df[(df["Species"] == "Yeast") & (df["alpha"] == 1.0)]["frac_closed"].mean()
    assert human < yeast, (
        f"Expected human frac_closed < yeast at α=1, "
        f"got human={human:.3f}, yeast={yeast:.3f}"
    )


def test_frac_closed_zero_at_alpha0(mock_full):
    """At α=0 (plain Vienna) the fraction closed must be exactly 0."""
    df = cofold_plot._add_frac_closed(mock_full)
    alpha0 = df[df["alpha"] == 0.0]["frac_closed"].to_numpy()
    assert np.allclose(alpha0, 0.0, atol=1e-10)


# ---------------------------------------------------------------------------
# API contracts
# ---------------------------------------------------------------------------

def test_gap_strip_panels_is_alias():
    assert cofold_plot.gap_strip_panels is cofold_plot.gap_heatmap_panels


def test_empty_dataframe_guard(tmp_path):
    assert cofold_plot.gap_closure_panels(pd.DataFrame(), tmp_path, "png") == []
    assert cofold_plot.gap_heatmap_panels(pd.DataFrame(), tmp_path, "png") == []
    assert cofold_plot.per_gene_landscape(pd.DataFrame(), tmp_path, "png") == []
