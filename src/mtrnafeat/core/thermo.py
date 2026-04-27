"""Thin ViennaRNA wrapper. All Vienna calls live here.

Functions return plain Python types, not Vienna objects, so the rest of the
package never imports `RNA` directly. Tests can mock this module wholesale.
"""
from __future__ import annotations

from typing import Optional

try:
    import RNA  # type: ignore
except ImportError as e:
    RNA = None  # type: ignore
    _IMPORT_ERROR = e
else:
    _IMPORT_ERROR = None


def _require_rna() -> None:
    if RNA is None:
        raise ImportError(
            "ViennaRNA Python bindings not installed. "
            "Install via `pip install viennarna` or `conda install -c bioconda viennarna`. "
            f"(original: {_IMPORT_ERROR})"
        )


def fold_mfe(seq: str, max_bp_span: Optional[int] = None) -> tuple[str, float]:
    """Return (dot-bracket, MFE in kcal/mol)."""
    _require_rna()
    if max_bp_span is None:
        struct, mfe = RNA.fold(seq)
        return struct, float(mfe)
    md = RNA.md()
    md.max_bp_span = int(max_bp_span)
    fc = RNA.fold_compound(seq, md)
    struct, mfe = fc.mfe()
    return struct, float(mfe)


def eval_structure(seq: str, structure: str, max_bp_span: Optional[int] = None) -> float:
    """Energy of a given dot-bracket structure under the standard NN model."""
    _require_rna()
    if max_bp_span is None:
        fc = RNA.fold_compound(seq)
    else:
        md = RNA.md()
        md.max_bp_span = int(max_bp_span)
        fc = RNA.fold_compound(seq, md)
    return float(fc.eval_structure(structure))


def ensemble_diversity(seq: str) -> float:
    """ViennaRNA mean base-pair distance across the Boltzmann ensemble.

    Uses `exp_params_rescale(mfe)` to avoid PF overflow on long sequences.
    """
    _require_rna()
    fc = RNA.fold_compound(seq)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    return float(fc.mean_bp_distance())


def centroid_structure(seq: str) -> tuple[str, float]:
    """Partition-function centroid structure and its distance from ensemble."""
    _require_rna()
    fc = RNA.fold_compound(seq)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    struct, dist = fc.centroid()
    return struct, float(dist)


def mea_structure(seq: str, gamma: float = 1.0) -> tuple[str, float]:
    """Maximum Expected Accuracy structure (gamma defaults to 1.0)."""
    _require_rna()
    fc = RNA.fold_compound(seq)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    struct, mea = fc.MEA(gamma)
    return struct, float(mea)


def bp_distance(struct_a: str, struct_b: str) -> int:
    """Hamming distance between two pair-tables (number of differing pairs)."""
    _require_rna()
    return int(RNA.bp_distance(struct_a, struct_b))
