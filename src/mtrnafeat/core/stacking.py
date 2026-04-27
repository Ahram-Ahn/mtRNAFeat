"""Turner-1999 dinucleotide stacking ΔG (kcal/mol) — sliding decomposition.

Lookup table is the canonical Watson-Crick / G-U wobble nearest-neighbor
stacking parameters from Mathews et al. 1999 (JMB 288:911-940), as
reproduced in ViennaRNA's `RNAfold -p` output.

Each entry maps a 5'-XY-3' / 3'-WZ-5' stacking step to its ΔG37 contribution.
We index by the top strand (5'-XY-3') and the bottom strand (3'-WZ-5'). When
asking for the stacking energy of two adjacent paired columns, we need the
identity of all four bases.

For convenience we expose a sliding-window function `stacking_track(seq,
struct)` that returns one ΔG value per stacking step where both adjacent
columns are paired.
"""
from __future__ import annotations

import numpy as np

from mtrnafeat.core.structure import pair_table

# Mathews 1999 (Table 4 / Vienna Turner 1999 set) values in kcal/mol.
# Key: ((top5, top3), (bot3, bot5)) where bottom is read 3'->5'.
# Equivalently: 5'-AB-3' paired with 3'-CD-5'.
# Source: Mathews et al. 1999 JMB; Vienna params/turner1999/rna_turner1999.par.
STACK37: dict[tuple[str, str, str, str], float] = {
    # AU/UA family
    ("A", "A", "U", "U"): -0.9,
    ("A", "U", "U", "A"): -1.1,
    ("U", "A", "A", "U"): -1.3,
    ("A", "C", "U", "G"): -2.2,
    ("A", "G", "U", "C"): -2.1,
    ("C", "A", "G", "U"): -2.1,
    ("C", "U", "G", "A"): -1.8,
    ("G", "A", "C", "U"): -2.4,
    ("G", "U", "C", "A"): -1.4,  # near-canonical, includes wobble end
    ("U", "C", "A", "G"): -2.3,
    ("U", "G", "A", "C"): -1.0,
    # GC/CG family
    ("C", "C", "G", "G"): -3.3,
    ("C", "G", "G", "C"): -2.4,
    ("G", "C", "C", "G"): -3.4,
    ("G", "G", "C", "C"): -3.3,
    # GU wobble family (selection, plus generic)
    ("G", "U", "U", "G"): -1.3,
    ("U", "G", "G", "U"): -0.6,
    ("G", "G", "U", "C"): -1.5,
    ("G", "C", "U", "G"): -2.5,
    ("C", "G", "G", "U"): -2.0,
    ("U", "G", "G", "C"): -1.5,
    ("G", "U", "C", "G"): -2.1,
    ("A", "G", "U", "U"): -1.4,
    ("G", "A", "U", "U"): -1.0,
    ("U", "A", "G", "U"): -1.0,
    ("U", "U", "A", "G"): -1.4,
    ("U", "U", "G", "A"): -0.9,
}


def stack_dG(t5: str, t3: str, b5: str, b3: str) -> float:
    """Return ΔG of the 5'-(t5)(t3)-3' / 3'-(b3)(b5)-5' step. 0.0 if missing.

    Note the bottom-strand direction: `b5` is the bottom base paired to t3
    (the 3'-most top base) on the bottom strand's 5' side; `b3` is the bottom
    base paired to t5.
    """
    key = (t5, t3, b3, b5)
    if key in STACK37:
        return STACK37[key]
    # Try reverse complement symmetry (canonical NN parameter symmetry):
    # 5'-XY-3'/3'-WZ-5' ≡ 5'-ZW-3'/3'-YX-5'
    rev = (b5, b3, t3, t5)
    return STACK37.get(rev, 0.0)


def stacking_track(seq: str, structure: str) -> np.ndarray:
    """Per-position stacking ΔG: at index i we report the ΔG of the step
    between paired columns i and i+1, when both are paired and to consecutive
    partners on the bottom strand. Otherwise 0.

    Returns an array of length len(seq)-1.
    """
    pt = pair_table(structure)
    n = len(seq)
    out = np.zeros(n - 1, dtype=float)
    for i in range(n - 1):
        if i in pt and (i + 1) in pt:
            j = pt[i]
            jp1 = pt[i + 1]
            if jp1 == j - 1:
                t5, t3 = seq[i], seq[i + 1]
                b5, b3 = seq[jp1], seq[j]
                out[i] = stack_dG(t5, t3, b5, b3)
    return out
