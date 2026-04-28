"""Subprocess wrapper for RNAstructure's ``Fold`` binary.

Defaults match Moran et al.: ``Fold -m 3 -md 350 [-dms <constraint>]``.
``fold()`` accepts an RNA sequence, an optional per-nucleotide DMS
reactivity vector (NaN = uninformative), and an optional ``max_bp_span``,
and returns a parsed dot-bracket + MFE.

Binary lookup (in priority order):
    1. ``RNASTRUCTURE_FOLD`` env var (full path to the binary)
    2. ``RNASTRUCTURE_BIN`` env var (directory containing ``Fold``)
    3. ``shutil.which("Fold")``

``DATAPATH`` must point at RNAstructure's ``data_tables/`` directory.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mtrnafeat.engines._common import (
    MissingEngineError,
    validate_rna_alphabet,
    write_dms_constraint_file,
)


@dataclass
class FoldResult:
    dot_bracket: str
    mfe: float


def _resolve_binary(name: str) -> str:
    override = os.environ.get(f"RNASTRUCTURE_{name.upper()}")
    if override:
        return override
    binroot = os.environ.get("RNASTRUCTURE_BIN")
    if binroot:
        cand = Path(binroot) / name
        if cand.exists():
            return str(cand)
    found = shutil.which(name)
    if found:
        return found
    raise MissingEngineError(
        f"RNAstructure '{name}' not on PATH. Install RNAstructure ≥6.4 and "
        f"export DATAPATH to its data_tables/ directory (or set "
        f"RNASTRUCTURE_BIN to the bin/ directory)."
    )


def _require_datapath() -> None:
    if "DATAPATH" not in os.environ:
        raise MissingEngineError(
            "RNAstructure DATAPATH env var is unset — point it at the "
            "data_tables/ directory shipped with RNAstructure."
        )


def _write_seq_file(sequence: str, path: Path, name: str = "seq") -> Path:
    """RNAstructure .seq format: ``;<comment>\\n<name>\\n<sequence>1`` (the
    trailing ``1`` is the canonical sentinel)."""
    with path.open("w") as fh:
        fh.write(";\n")
        fh.write(f"{name}\n")
        fh.write(f"{sequence}1\n")
    return path


def _ct_to_dot_bracket(ct_path: Path) -> tuple[str, float]:
    """Parse a CT file and return (dot-bracket, energy_kcal_per_mol).

    Only the first structure record is read (Fold -m 1 default). The header is
    ``<N>  ENERGY = <e>  <name>``; rows are ``<i> <base> <i-1> <i+1> <pair> <i>``.
    """
    with ct_path.open() as fh:
        header = fh.readline()
        m = re.search(r"ENERGY\s*=\s*(-?\d+\.\d+)", header)
        energy = float(m.group(1)) if m else float("nan")
        try:
            n = int(header.split()[0])
        except (ValueError, IndexError) as exc:
            raise RuntimeError(f"Malformed CT header: {header!r}") from exc
        pair = [0] * (n + 1)
        for line in fh:
            parts = line.split()
            if len(parts) < 6:
                continue
            i = int(parts[0])
            j = int(parts[4])
            if i > n:  # next structure block — stop
                break
            pair[i] = j
        out = ["."] * n
        for i in range(1, n + 1):
            j = pair[i]
            if j and j > i:
                out[i - 1] = "("
                out[j - 1] = ")"
        return "".join(out), energy


def _run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess:
    proc = subprocess.run(
        cmd, cwd=str(cwd), capture_output=True, text=True, check=False
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"{cmd[0]} failed (rc={proc.returncode}): {proc.stderr.strip()}"
        )
    return proc


def fold(
    seq: str,
    dms: np.ndarray | None = None,
    max_bp_span: int = 350,
) -> FoldResult:
    """RNAstructure ``Fold`` — MFE-style DMS-guided structure.

    ``-m 3`` keeps the top 3 suboptimal structures internally; we read only
    the first record. ``-md`` caps base-pair distance to ``max_bp_span``.
    Pass ``dms`` to add a ``-dms <constraint>`` flag; otherwise the fold is
    purely thermodynamic.
    """
    validate_rna_alphabet(seq)
    bin_path = _resolve_binary("Fold")
    _require_datapath()
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        seq_path = _write_seq_file(seq, td_path / "in.seq")
        ct_path = td_path / "out.ct"
        cmd = [
            bin_path, str(seq_path), str(ct_path),
            "-m", "3",
            "-md", str(int(max_bp_span)),
        ]
        if dms is not None:
            dms_path = write_dms_constraint_file(dms, seq, td_path / "in.dms")
            cmd.extend(["-dms", str(dms_path)])
        _run(cmd, cwd=td_path)
        dot, energy = _ct_to_dot_bracket(ct_path)
    return FoldResult(dot_bracket=dot, mfe=energy)
