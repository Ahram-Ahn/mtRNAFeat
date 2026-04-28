"""Canonical names, region labels, and palette. Single source of truth."""
from __future__ import annotations

# Resolves the .db file naming (ATP86, ND4L4) to the codon_data annotation
# naming (ATP8/ATP6, ND4L/ND4) and back. All readers normalize on entry.
GENE_ALIASES: dict[str, str] = {
    "ATP86": "ATP8/ATP6",
    "ATP8/ATP6": "ATP8/ATP6",
    "ATP8_ATP6": "ATP8/ATP6",
    "ND4L4": "ND4L/ND4",
    "ND4L/ND4": "ND4L/ND4",
    "ND4L_ND4": "ND4L/ND4",
}

# The full mt-mRNA gene set we expect across both species (.db files).
EXPECTED_GENES = (
    "COX1", "COX2", "COX3",
    "CYTB", "COB",
    "ATP86", "ATP9",
    "ND1", "ND2", "ND3", "ND4L4", "ND5", "ND6",
    "VAR1",
)

REGION_LABELS = ("5'UTR", "CDS", "3'UTR/tail")

# Publication palette — used by every viz module via viz.style.apply().
PALETTE = {
    "Human": "#D62728",      # red
    "Yeast": "#FF7F0E",      # orange
    "Sim Human": "#F4A6A2",
    "Sim Yeast": "#FFD7A8",
    "DMS": "#1F77B4",        # blue
    "Vienna": "#D62728",
    "Vienna_span": "#2CA02C",
    "RNAstructure": "#8C564B",   # brown — distinct from the Human/Vienna red

    "Diversity": "#000000",
    "Cotrans": "#9467BD",
    "Difference": "#7F7F7F",
    "UTR": "#BDBDBD",
    "CDS": "#4DAF4A",
}

# Fixed CSV float format so determinism tests can SHA-compare.
CSV_FLOAT_FORMAT = "%.10g"

# What the analysis modules emit in their KIND header so downstream code
# never confuses three different kinds of "co-transcriptional" computation.
KIND_PROJECTED = "projected_realization"
KIND_PREFIX_DIVERSITY = "ensemble_prefix_diversity"
KIND_TRUE_KINETIC = "true_kinetic"


def canonical_gene(name: str) -> str:
    """Return the canonical (annotation-table) form of a gene name."""
    n = name.strip()
    return GENE_ALIASES.get(n, n)


def db_gene(name: str) -> str:
    """Return the .db file form of a gene name (the inverse direction)."""
    canon = canonical_gene(name)
    inverse = {"ATP8/ATP6": "ATP86", "ND4L/ND4": "ND4L4"}
    return inverse.get(canon, canon)


def file_safe_gene(name: str) -> str:
    """Canonical gene name with path separators replaced for safe filenames."""
    return canonical_gene(name).replace("/", "_")
