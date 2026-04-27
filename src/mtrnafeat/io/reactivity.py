"""Loaders for per-nucleotide DMS/SHAPE reactivity profiles.

Stub: the user does not yet have raw reactivity data — only DMS-guided
dot-bracket structures. This module is reserved for when ShapeMapper2
`_profile.txt` or RNA Framework XML output becomes available, so that
ViennaRNA soft constraints can be activated downstream.
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd


def load_shapemapper_profile(path: str | Path) -> pd.DataFrame:
    """Load a ShapeMapper2 _profile.txt file.

    Returns a DataFrame indexed by 1-based nucleotide position with at
    least 'reactivity' and 'sequence' columns. Raises NotImplementedError
    until the user provides a concrete sample profile to test against.
    """
    raise NotImplementedError(
        "Reactivity ingestion not yet wired. Provide a ShapeMapper2 _profile.txt "
        "fixture in test_data/ and unblock this loader."
    )
