"""Determinism: re-running with the same seed reproduces byte-identical CSVs.

Skipped when ViennaRNA is not installed.
"""
from __future__ import annotations

import hashlib
import subprocess
import sys

from tests.conftest import needs_rna


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


@needs_rna
def test_landscape_deterministic(mini_config_path, tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    for outdir in (a, b):
        cmd = [sys.executable, "-m", "mtrnafeat.cli", "landscape",
               "--config", str(mini_config_path), "--outdir", str(outdir)]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        assert proc.returncode == 0, proc.stderr
    csv_a = a / "landscape" / "specific_conditions.csv"
    csv_b = b / "landscape" / "specific_conditions.csv"
    assert _sha256(csv_a) == _sha256(csv_b)
