"""Shared helpers for engine wrappers: alphabet validation, missing-engine
errors, and the RNAstructure DMS constraint-file format.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

_RNA_ALPHABET = frozenset("ACGU")


class MissingEngineError(RuntimeError):
    """Raised when a required external folding engine is unavailable.

    Distinguishes "the user forgot to install RNAstructure" from "the call
    crashed mid-fold" so callers can decide whether to skip or abort.
    """


def validate_rna_alphabet(seq: str) -> None:
    """Reject anything that isn't ACGU (case-sensitive). T must be converted
    to U upstream — the engine wrappers refuse to do it silently.
    """
    if not seq:
        raise ValueError("Empty sequence")
    bad = [(i + 1, ch) for i, ch in enumerate(seq) if ch not in _RNA_ALPHABET]
    if not bad:
        return
    first = bad[0]
    if first[1].upper() == "T":
        raise ValueError(
            f"T detected at position {first[0]} — convert to U before folding "
            f"(engines refuse silent T→U conversion)"
        )
    if first[1] in "acgu":
        raise ValueError(
            f"Lowercase nucleotide '{first[1]}' at position {first[0]} — "
            f"upper-case the sequence before folding"
        )
    raise ValueError(
        f"Non-RNA character '{first[1]}' at position {first[0]} "
        f"(allowed: A,C,G,U)"
    )


def dms_to_minus999(dms: np.ndarray, sequence: str) -> np.ndarray:
    """Translate a per-nucleotide reactivity array to RNAstructure's wire format.

    Input: float array of length len(sequence). NaN means "no information"
    (low coverage or chemically masked). Output: same length, with -999 at
    every G/U position (DMS doesn't react with G/U) and at every NaN position.
    """
    if len(dms) != len(sequence):
        raise ValueError(
            f"Reactivity length {len(dms)} != sequence length {len(sequence)}"
        )
    out = np.asarray(dms, dtype=float).copy()
    for i, base in enumerate(sequence):
        if base in ("G", "U") or not np.isfinite(out[i]):
            out[i] = -999.0
    return out


def write_dms_constraint_file(
    dms: np.ndarray, sequence: str, path: Path
) -> Path:
    """RNAstructure DMS/SHAPE constraint file: tab-separated `index value`
    pairs, no header, 1-based positions, -999 for any uninformative position.
    """
    arr = dms_to_minus999(dms, sequence)
    path = Path(path)
    with path.open("w") as fh:
        for i, v in enumerate(arr, start=1):
            fh.write(f"{i}\t{v:.6f}\n")
    return path
