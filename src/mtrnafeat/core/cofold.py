"""CoFold-style co-transcriptional folding (pure Python, ViennaRNA soft constraints).

Inspired by the co-transcriptional folding approach of Proctor & Meyer
2013, *NAR* Vol 41 Iss 9, e102
(https://academic.oup.com/nar/article/41/9/e102/2409154 ;
DOI 10.1093/nar/gkt174). Implemented here as a kcal/mol soft long-range
penalty:

    f(d) = alpha * (1 - exp(-d / tau))

applied as an additional energy on each candidate base pair (i, j) with
sequence distance d = j - i. Long-range pairs are progressively
penalized, modelling the kinetic disadvantage of pairing partners that
were transcribed far apart in time.

NOTE on parameter values: the published CoFold paper uses a dimensionless
multiplicative scaling on candidate-pair "reachability" with very
different numerical defaults. The mtrnafeat defaults alpha = 0.5
kcal/mol and tau = 640 nt are this implementation's own choices: alpha
is in kcal/mol (not dimensionless), and tau = 640 corresponds to a
~50 nt/s transcription rate and ~12.8 s pairing window under that
parameterization. Treat sweep outputs as exploratory rather than as a
literal reproduction of the published CoFold method.

We implement this purely via ViennaRNA's `fc.sc_add_bp` interface — no
external binary, no patched libRNA. At alpha=0 the result is exactly
identical to the standard MFE.

This module is the default thermodynamic engine for window-level and
TIS-level analyses. Plain Vienna MFE is still emitted side-by-side so
results can be compared directly.
"""
from __future__ import annotations

import numpy as np

from mtrnafeat.core.thermo import _require_rna


def _penalty_kcal(distance: int, alpha: float, tau: float) -> float:
    """CoFold soft penalty for a pair with sequence distance `distance`, in kcal/mol."""
    if distance <= 0:
        return 0.0
    return float(alpha * (1.0 - np.exp(-distance / tau)))


def cofold_dG(seq: str, alpha: float = 0.5, tau: float = 640.0,
              max_bp_span: int | None = None,
              min_pair_distance: int = 4) -> tuple[str, float]:
    """Fold `seq` under the CoFold soft long-range penalty.

    Parameters
    ----------
    seq : RNA sequence (A/C/G/U/T; T normalized in caller).
    alpha : penalty strength in kcal/mol; CoFold paper uses 0.5.
    tau : decay constant in nt; CoFold paper uses 640.
    max_bp_span : optional hard cutoff (no pairs allowed beyond this distance).
        Useful in combination with the soft penalty for the window-scan use case.
    min_pair_distance : minimum sequence distance for a valid pair (RNA min
        hairpin loop = 3, so this defaults to 4 = i,j with j-i >= 4).

    Returns
    -------
    (dot_bracket, mfe_kcal): the optimal structure under the modified model
        and its energy. Note the energy returned by ViennaRNA already includes
        the soft-constraint penalties on whichever pairs are in the optimum.
    """
    _require_rna()
    import RNA  # noqa: WPS433
    n = len(seq)
    if n < min_pair_distance + 1:
        return "." * n, 0.0

    if max_bp_span is None:
        fc = RNA.fold_compound(seq)
    else:
        md = RNA.md()
        md.max_bp_span = int(max_bp_span)
        fc = RNA.fold_compound(seq, md)

    if alpha != 0.0:
        for i in range(1, n):
            for j in range(i + min_pair_distance, n + 1):
                d = j - i
                if max_bp_span is not None and d > int(max_bp_span):
                    break  # j-loop is monotone, can break early
                pen = _penalty_kcal(d, alpha, tau)
                if pen > 0.0:
                    fc.sc_add_bp(i, j, pen)

    structure, mfe = fc.mfe()
    return structure, float(mfe)


def cofold_eval(seq: str, structure: str, alpha: float = 0.5, tau: float = 640.0,
                 max_bp_span: int | None = None) -> float:
    """Energy of `structure` under the CoFold-modified scoring (kcal/mol).

    Computes Vienna's plain eval_structure plus the sum of CoFold penalties
    over all base pairs in `structure`. This avoids re-attaching soft
    constraints just to evaluate, which is much faster for batch scoring.
    """
    _require_rna()
    import RNA  # noqa: WPS433
    if max_bp_span is None:
        fc = RNA.fold_compound(seq)
    else:
        md = RNA.md()
        md.max_bp_span = int(max_bp_span)
        fc = RNA.fold_compound(seq, md)
    base = float(fc.eval_structure(structure))
    if alpha == 0.0:
        return base
    pairs = _extract_pairs(structure)
    penalty_kcal = 0.0
    for i, j in pairs:
        d = (j - i)
        penalty_kcal += alpha * (1.0 - np.exp(-d / tau))
    return base + penalty_kcal


def _extract_pairs(structure: str) -> list[tuple[int, int]]:
    stack: list[int] = []
    pairs: list[tuple[int, int]] = []
    for k, ch in enumerate(structure, start=1):
        if ch == "(":
            stack.append(k)
        elif ch == ")":
            if stack:
                i = stack.pop()
                pairs.append((i, k))
    return pairs
