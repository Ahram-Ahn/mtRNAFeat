"""Per-transcript UTR/CDS annotations. Replaces legacy/codon_data.py annotation_df.

Gene names are stored in canonical form (ATP8/ATP6, ND4L/ND4) and aliased to
their .db form (ATP86, ND4L4) by `mtrnafeat.constants.canonical_gene`.
"""
from __future__ import annotations

import pandas as pd

from mtrnafeat.constants import canonical_gene

_YEAST = pd.DataFrame({
    "transcript": ["COX1", "ATP8/ATP6", "COB", "ATP9", "VAR1", "COX2", "COX3"],
    "l_tr":   [2147, 2065, 2223,  960, 1952,  886, 1528],
    "l_utr5": [ 460,  364,  954,  630,  162,   54,  604],
    "l_utr3": [  82,  100,  111,   99,  593,   76,  114],
})
_YEAST["l_cds"] = _YEAST["l_tr"] - _YEAST["l_utr5"] - _YEAST["l_utr3"]

_HUMAN = pd.DataFrame({
    "transcript": ["ND1", "ND2", "COX1", "COX2", "ATP8/ATP6", "COX3",
                   "ND3", "ND4L/ND4", "ND5", "ND6", "CYTB"],
    "l_tr":   [958, 1042, 1617,  708,  843,  784,  346, 1668, 2379,  538, 1141],
    "l_utr5": [  2,    0,    3,    0,    1,    0,    0,    0,    0,    0,    0],
    "l_utr3": [  0,    0,   72,   24,    0,    0,    0,    0,  567,   13,    0],
})
_HUMAN["l_cds"] = _HUMAN["l_tr"] - _HUMAN["l_utr5"] - _HUMAN["l_utr3"]


def annotation_df(species: str) -> pd.DataFrame:
    """Return the annotation DataFrame for a species (Human/Yeast)."""
    s = species.lower()
    if s == "human":
        return _HUMAN.copy()
    if s == "yeast":
        return _YEAST.copy()
    raise ValueError(f"Unknown species: {species}")


def annotation_for(species: str, gene: str) -> dict[str, int]:
    """Return {l_tr, l_utr5, l_cds, l_utr3} for one transcript."""
    df = annotation_df(species)
    target = canonical_gene(gene)
    row = df[df["transcript"] == target]
    if row.empty:
        raise KeyError(f"{species}: no annotation for {gene} (canonical {target})")
    r = row.iloc[0]
    return {
        "l_tr": int(r["l_tr"]),
        "l_utr5": int(r["l_utr5"]),
        "l_cds": int(r["l_cds"]),
        "l_utr3": int(r["l_utr3"]),
    }


def classify_region(l_utr5: int, l_cds: int, pos1: int) -> str:
    """1-based position → region label."""
    if pos1 <= l_utr5:
        return "5'UTR"
    if pos1 <= l_utr5 + l_cds:
        return "CDS"
    return "3'UTR/tail"
