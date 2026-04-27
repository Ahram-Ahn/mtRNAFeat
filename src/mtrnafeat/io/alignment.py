"""Parser for the PAL2NAL aa-dna alignment of yeast/human COX1.

PAL2NAL output is verbose; for our purposes we only need the codon-aligned
yeast and human DNA. Translations are recomputed from codons via the
mitochondrial codon tables in `mtrnafeat.io.codons`.

We extract codon strings from every line that starts with `S.cerevisiae` or
`H.sapiens` and pair them up in alignment order.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from mtrnafeat.io.codons import HUMAN_MT_CODON_TABLE, YEAST_MT_CODON_TABLE


@dataclass
class CodonAlignment:
    yeast_codons: list[str]
    human_codons: list[str]
    yeast_aa: list[str]
    human_aa: list[str]

    def __post_init__(self) -> None:
        n = len(self.yeast_codons)
        if not (len(self.human_codons) == len(self.yeast_aa) == len(self.human_aa) == n):
            raise ValueError("Alignment columns out of sync")

    @property
    def n_columns(self) -> int:
        return len(self.yeast_codons)


def _translate(codon: str, table: dict[str, str]) -> str:
    if "-" in codon:
        return "-"
    return table.get(codon.upper(), "?")


def parse_pal2nal(path: str | Path) -> CodonAlignment:
    yeast_codons: list[str] = []
    human_codons: list[str] = []
    with open(path) as fh:
        for ln in fh:
            stripped = ln.rstrip("\n")
            if stripped.lstrip().startswith("S.cerevisiae"):
                yeast_codons.extend(_codons_from(stripped, "S.cerevisiae"))
            elif stripped.lstrip().startswith("H.sapiens"):
                human_codons.extend(_codons_from(stripped, "H.sapiens"))
    n = min(len(yeast_codons), len(human_codons))
    yeast_codons = yeast_codons[:n]
    human_codons = human_codons[:n]
    yeast_aa = [_translate(c, YEAST_MT_CODON_TABLE) for c in yeast_codons]
    human_aa = [_translate(c, HUMAN_MT_CODON_TABLE) for c in human_codons]
    return CodonAlignment(yeast_codons=yeast_codons, human_codons=human_codons,
                          yeast_aa=yeast_aa, human_aa=human_aa)


def _codons_from(line: str, label: str) -> list[str]:
    rest = line.split(label, 1)[1]
    return [tok for tok in rest.split() if len(tok) == 3]


def codon_position_changes(codon_y: str, codon_h: str) -> list[int]:
    if "-" in codon_y or "-" in codon_h:
        return []
    return [i + 1 for i in range(3) if codon_y[i] != codon_h[i]]
