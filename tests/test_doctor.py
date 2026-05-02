"""Smoke tests for `mtrnafeat doctor`."""
from __future__ import annotations

import json
import subprocess
import sys


def _run(*extra: str) -> subprocess.CompletedProcess:
    cmd = [sys.executable, "-m", "mtrnafeat.cli", "doctor", *extra]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=30)


def test_doctor_runs_and_prints_summary():
    """Doctor must execute end-to-end and print the standard label rows.

    Exit code can be 0 or 1 depending on the test environment (e.g. RNA
    binding installed or not). What we assert is that it completes
    without raising and prints the labels we promised in STAGES.md, so
    that script-driven users can grep for them.
    """
    proc = _run()
    assert proc.returncode in (0, 1), proc.stderr
    for label in (
        "Python >=3.11",
        "mtrnafeat import",
        "ViennaRNA import",
        "RNAplfold availability",
        "RNAstructure Fold",
        "DrTransformer (kinetic)",
    ):
        assert label in proc.stdout, f"missing label {label!r} in stdout:\n{proc.stdout}"


def test_doctor_json_output(tmp_path):
    out_json = tmp_path / "doctor.json"
    proc = _run("--", "--json", str(out_json))
    assert proc.returncode in (0, 1), proc.stderr
    assert out_json.exists()
    data = json.loads(out_json.read_text())
    assert isinstance(data, list) and data, "JSON output must be a non-empty list"
    sample = data[0]
    assert set(sample.keys()) == {"level", "item", "message"}
    assert sample["level"] in {"OK", "WARN", "ERROR", "OPTIONAL"}
