"""Tests for run-all orchestration behavior."""
from __future__ import annotations

from types import SimpleNamespace

from mtrnafeat.commands import pipeline
from mtrnafeat.config import Config


def test_sequential_returns_nonzero_on_stage_failure(monkeypatch, tmp_path):
    cfg = Config(outdir=tmp_path, db_files={"Sample": "unused.db"})

    def fake_import_module(name: str):
        if name.endswith(".bad"):
            return SimpleNamespace(run=lambda _cfg, _args: 1)
        return SimpleNamespace(run=lambda _cfg, _args: 0)

    monkeypatch.setattr(pipeline.importlib, "import_module", fake_import_module)
    rc = pipeline._sequential(cfg, skip=set(), stages=("good", "bad"))
    assert rc == 1


def test_parallel_child_stages_force_single_worker(monkeypatch, tmp_path):
    cfg = Config(outdir=tmp_path, db_files={"Sample": "unused.db"}, seed=7)
    calls = []

    class DummyProc:
        returncode = 0

    def fake_run(cmd, capture_output, env):
        calls.append((cmd, env))
        return DummyProc()

    monkeypatch.setattr(pipeline.subprocess, "run", fake_run)
    rc = pipeline._parallel(cfg, skip=set(), config_path=None, stages=("stats",))
    assert rc == 0
    assert len(calls) == 1
    assert calls[0][1]["MTRNAFEAT_N_WORKERS"] == "1"


def test_parallel_returns_nonzero_on_stage_failure(monkeypatch, tmp_path):
    cfg = Config(outdir=tmp_path, db_files={"Sample": "unused.db"}, seed=7)

    class DummyProc:
        returncode = 2

    monkeypatch.setattr(pipeline.subprocess, "run", lambda cmd, capture_output, env: DummyProc())
    rc = pipeline._parallel(cfg, skip=set(), config_path=None, stages=("stats",))
    assert rc == 1
