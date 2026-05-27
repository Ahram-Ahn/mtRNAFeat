"""Substitution-thermo bookkeeping tests that do not require ViennaRNA."""
from __future__ import annotations

import os
import subprocess
import sys

import pandas as pd
import pytest

from mtrnafeat.analysis.substitution_thermo import _stable_species_seed_offset, _summarize


def test_species_seed_offset_is_not_python_hash_dependent():
    code = (
        "from mtrnafeat.analysis.substitution_thermo import _stable_species_seed_offset; "
        "print(_stable_species_seed_offset('Human'), _stable_species_seed_offset('Yeast'))"
    )
    vals = []
    for hash_seed in ("1", "2"):
        env = {**os.environ, "PYTHONHASHSEED": hash_seed}
        proc = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
        vals.append(proc.stdout.strip())
    assert vals[0] == vals[1] == "101 211"
    assert _stable_species_seed_offset("Human") == 101


def test_summary_reports_chunk_and_full_cds_energy_scale():
    dist = pd.DataFrame(
        [
            {
                "Species": "Human",
                "Gene": "ND2",
                "Pool": "WildType_MFE",
                "MFE_kcal_per_mol": -30.0,
                "Length_nt": 300,
            },
            {
                "Species": "Human",
                "Gene": "ND2",
                "Pool": "WildType_DMS_Eval",
                "MFE_kcal_per_mol": -3.0,
                "Length_nt": 300,
            },
            {
                "Species": "Human",
                "Gene": "ND2",
                "Pool": "WildType_DMS_Eval_FullLength",
                "MFE_kcal_per_mol": -60.0,
                "Length_nt": 1000,
            },
            *[
                {
                    "Species": "Human",
                    "Gene": "ND2",
                    "Pool": "synonymous",
                    "MFE_kcal_per_mol": mfe,
                    "Length_nt": 300,
                }
                for mfe in (-40.0, -35.0, -30.0)
            ],
        ]
    )

    summary = _summarize(dist)
    row = summary[summary["Pool"] == "synonymous"].iloc[0]

    assert row["Length_nt"] == 300
    assert row["DMS_Full_CDS_Length_nt"] == 1000
    assert row["WT_DMS_Eval_per_nt"] == pytest.approx(-0.01)
    assert row["WT_DMS_Eval_FullLength_per_nt"] == pytest.approx(-0.06)
    assert "codon-aligned CDS prefix" in row["DMS_dG_Source"]
