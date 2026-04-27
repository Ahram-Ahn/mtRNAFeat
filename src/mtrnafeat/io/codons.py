"""Codon tables. Standard, yeast-mt, and human-mt."""
from __future__ import annotations

# Lifted verbatim from legacy/codon_data.py.

STANDARD_CODON_TABLE: dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

YEAST_MT_CODON_TABLE: dict[str, str] = {
    **STANDARD_CODON_TABLE,
    'TGA': 'W',
    'CTT': 'T', 'CTC': 'T', 'CTA': 'T', 'CTG': 'T',
    'ATA': 'M',
}

HUMAN_MT_CODON_TABLE: dict[str, str] = {
    **STANDARD_CODON_TABLE,
    'TGA': 'W',
    'ATA': 'M',
    'AGA': '*', 'AGG': '*',
}


def codon_table_for(species: str) -> dict[str, str]:
    s = species.lower()
    if s == "human":
        return HUMAN_MT_CODON_TABLE
    if s == "yeast":
        return YEAST_MT_CODON_TABLE
    if s == "standard":
        return STANDARD_CODON_TABLE
    raise ValueError(f"Unknown species: {species}")
