"""Parser for the 3-line .db format produced by the DMS-MaPseq pipeline.

Each record:
    >GENE: -X kcal/mol
    SEQUENCE
    DOTBRACKET
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from mtrnafeat.constants import canonical_gene


@dataclass(frozen=True)
class DbRecord:
    gene: str         # canonical name
    raw_gene: str     # name as written in the .db file
    mfe: float
    sequence: str     # always upper-case RNA (T -> U)
    structure: str

    def __post_init__(self) -> None:
        if len(self.sequence) != len(self.structure):
            raise ValueError(
                f"{self.raw_gene}: sequence length {len(self.sequence)} != structure length {len(self.structure)}"
            )
        _validate_brackets(self.structure, self.raw_gene)


def _validate_brackets(struct: str, gene: str) -> None:
    depth = 0
    for ch in struct:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                raise ValueError(f"{gene}: unbalanced ')' in structure")
    if depth != 0:
        raise ValueError(f"{gene}: {depth} unclosed '(' in structure")


def _normalize_seq(seq: str) -> str:
    return seq.strip().upper().replace("T", "U")


def _parse_mfe(header: str) -> float:
    parts = header.split(":", 1)
    if len(parts) != 2:
        return 0.0
    try:
        return float(parts[1].replace("kcal/mol", "").strip())
    except ValueError:
        return 0.0


def parse_db(path: str | Path) -> list[DbRecord]:
    """Read a .db file and return one DbRecord per gene. Order preserved."""
    path = Path(path)
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    records: list[DbRecord] = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.startswith(">"):
            header = line[1:]
            raw_name = header.split(":", 1)[0].strip()
            mfe = _parse_mfe(header)
            if i + 2 >= len(lines):
                raise ValueError(f"Incomplete record for {raw_name} in {path}")
            seq = _normalize_seq(lines[i + 1])
            struct = lines[i + 2].strip()
            records.append(DbRecord(
                gene=canonical_gene(raw_name),
                raw_gene=raw_name,
                mfe=mfe,
                sequence=seq,
                structure=struct,
            ))
            i += 3
        else:
            i += 1
    return records


def get_record(path: str | Path, gene: str) -> DbRecord:
    """Return one record by canonical or .db-form gene name."""
    target = canonical_gene(gene)
    for rec in parse_db(path):
        if rec.gene == target:
            return rec
    raise KeyError(f"Gene {gene} not found in {path}")


def list_genes(path: str | Path) -> list[str]:
    """Return canonical gene names present in a .db file, in file order."""
    return [rec.gene for rec in parse_db(path)]
