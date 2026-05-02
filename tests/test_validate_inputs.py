"""Tests for `mtrnafeat validate-inputs` and the underlying validation primitives."""
from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from mtrnafeat.validation import (
    check_db_file,
    has_error,
    run_all_checks,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _run_cli(*extra: str) -> subprocess.CompletedProcess:
    cmd = [sys.executable, "-m", "mtrnafeat.cli", "validate-inputs", *extra]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=30)


# ---------------------------------------------------------------------------
# Direct primitive tests
# ---------------------------------------------------------------------------

def test_check_db_length_mismatch_is_error():
    issues = check_db_file(FIXTURES / "bad_length.db")
    assert has_error(issues), [i for i in issues]
    assert any("length" in i.message.lower() for i in issues)


def test_check_db_unbalanced_brackets_is_error():
    issues = check_db_file(FIXTURES / "bad_brackets.db")
    assert has_error(issues), [i for i in issues]
    assert any("unbalanced" in i.message.lower() or "unclosed" in i.message.lower()
               for i in issues)


def test_check_db_missing_file_is_error(tmp_path):
    issues = check_db_file(tmp_path / "does_not_exist.db")
    assert has_error(issues)


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------

def test_validate_inputs_valid_config(mini_config_path):
    """Mini fixture is a known-good config; expect exit 0."""
    proc = _run_cli("--config", str(mini_config_path))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "OK" in proc.stdout


def test_validate_inputs_no_config_fails():
    proc = _run_cli()
    assert proc.returncode == 1
    assert "requires --config" in proc.stderr


def test_validate_inputs_json_output(mini_config_path, tmp_path):
    out = tmp_path / "v.json"
    proc = _run_cli("--config", str(mini_config_path), "--", "--json", str(out))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    data = json.loads(out.read_text())
    assert isinstance(data, list) and data
    assert all(set(row) == {"level", "item", "message"} for row in data)


def test_validate_inputs_bad_db(tmp_path):
    """Point a config at a bad .db file and expect non-zero exit."""
    bad_dir = tmp_path / "data"
    bad_dir.mkdir()
    (bad_dir / "human.db").write_text((FIXTURES / "bad_length.db").read_text())
    (bad_dir / "yeast.db").write_text((FIXTURES / "bad_brackets.db").read_text())
    cfg_path = tmp_path / "bad.yaml"
    cfg_path.write_text(textwrap.dedent(f"""\
        seed: 42
        data_dir: {bad_dir}
        outdir: {tmp_path / 'runs'}
        db_files:
          Human: human.db
          Yeast: yeast.db
        target_genes: [COX1]
        fold_engine: vienna
    """))

    proc = _run_cli("--config", str(cfg_path))
    assert proc.returncode == 1
    assert "ERROR" in proc.stdout


def test_validate_inputs_missing_optional_alignment(tmp_path, mini_config_path):
    """Missing alignment file must be a WARN, not an ERROR."""
    issues = run_all_checks(mini_config_path)
    # mini_config_path doesn't ship an alignment; we expect a WARN row but
    # no ERROR rows from the alignment check itself.
    alignment_rows = [i for i in issues if i.item == "alignment"]
    assert alignment_rows, "alignment row missing from validation output"
    assert all(i.level != "ERROR" for i in alignment_rows)


@pytest.mark.parametrize("fixture_name,expect_error", [
    ("bad_length.db", True),
    ("bad_brackets.db", True),
])
def test_check_db_fixtures_parametric(fixture_name, expect_error):
    issues = check_db_file(FIXTURES / fixture_name)
    assert has_error(issues) is expect_error
