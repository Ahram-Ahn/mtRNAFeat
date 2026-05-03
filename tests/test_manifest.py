"""Tests for the run-manifest writer."""
from __future__ import annotations

import json
from pathlib import Path

from mtrnafeat.core.manifest import write_run_manifest


def test_write_run_manifest_includes_required_keys(tmp_path: Path):
    path = write_run_manifest(
        tmp_path,
        command="stats",
        argv=["stats", "--config", "configs/all.yaml"],
        config_path=Path("configs/all.yaml"),
        seed=42,
    )
    assert path.exists()
    data = json.loads(path.read_text())
    for key in (
        "mtrnafeat_version", "command", "argv", "config_path", "seed",
        "viennarna_version", "rnastructure_version", "git_commit",
        "python_version", "platform", "timestamp_utc",
    ):
        assert key in data, f"missing key {key!r} in manifest"
    assert data["command"] == "stats"
    assert data["seed"] == 42
    assert data["argv"] == ["stats", "--config", "configs/all.yaml"]


def test_write_run_manifest_overwrites_previous(tmp_path: Path):
    write_run_manifest(
        tmp_path, command="stats", argv=["stats"],
        config_path=None, seed=1,
    )
    p = write_run_manifest(
        tmp_path, command="landscape", argv=["landscape"],
        config_path=None, seed=2,
    )
    data = json.loads(p.read_text())
    assert data["command"] == "landscape"
    assert data["seed"] == 2
