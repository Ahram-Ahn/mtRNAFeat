"""Kinetic-stage robustness tests."""
from __future__ import annotations

from mtrnafeat.commands import kinetic as kinetic_command
from mtrnafeat.config import Config


def test_kinetic_runtime_failure_returns_optional_dependency_code(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(kinetic_command.kinetic, "has_drtransformer", lambda: True)

    def fail_run(_cfg, _genes):
        raise RuntimeError("DrTransformer failed: missing RNA")

    monkeypatch.setattr(kinetic_command.kinetic, "run_kinetic_for_genes", fail_run)
    cfg = Config(outdir=tmp_path)

    rc = kinetic_command.run(cfg, [])

    assert rc == 2
    assert "kinetic unavailable or failed" in capsys.readouterr().out
