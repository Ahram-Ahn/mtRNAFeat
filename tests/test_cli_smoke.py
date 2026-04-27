"""Smoke test: every subcommand runs cleanly on the mini fixture in <30s.

Skipped when ViennaRNA is not installed.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from tests.conftest import needs_rna


SUBCOMMANDS = [
    ("stats", []),
    ("landscape", []),
    ("features", []),
    ("window", []),
    ("significance", []),
    ("tis", []),
    ("compare", []),
    ("substitution", ["--n", "5"]),
    ("cofold", []),
    ("gene-panel", []),
]


@needs_rna
@pytest.mark.parametrize("name,extra", SUBCOMMANDS)
def test_subcommand_smoke(name, extra, mini_config_path, tmp_path):
    cmd = [sys.executable, "-m", "mtrnafeat.cli", name, "--config", str(mini_config_path),
           "--outdir", str(tmp_path)]
    if extra:
        cmd += ["--", *extra]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, f"{name} failed:\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
