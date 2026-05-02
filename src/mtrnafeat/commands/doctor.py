"""`mtrnafeat doctor` — environment diagnostics.

Quickly tells the user whether the package and its external dependencies
are usable. Required components produce ERROR; optional/exploratory
components produce WARN or OPTIONAL. Exit status:

  - 0 if no ERRORs (warnings are tolerated).
  - 1 if any required component fails.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path

from mtrnafeat import __version__
from mtrnafeat.config import Config
from mtrnafeat.validation import (
    ValidationIssue,
    format_table,
    has_error,
    issues_to_rows,
    issues_writable,
    serialize_issues,
)


def _check_python() -> ValidationIssue:
    v = sys.version_info
    label = f"{v.major}.{v.minor}.{v.micro}"
    if (v.major, v.minor) < (3, 11):
        return ValidationIssue("ERROR", "Python >=3.11", label)
    return ValidationIssue("OK", "Python >=3.11", label)


def _check_self() -> ValidationIssue:
    return ValidationIssue("OK", "mtrnafeat import", f"version {__version__}")


def _check_module(modname: str, label: str) -> ValidationIssue:
    spec = importlib.util.find_spec(modname)
    if spec is None:
        return ValidationIssue("ERROR", label, f"module {modname!r} not importable")
    try:
        mod = importlib.import_module(modname)
    except Exception as exc:  # pragma: no cover — module import edge cases
        return ValidationIssue("ERROR", label, f"import failed: {type(exc).__name__}: {exc}")
    version = getattr(mod, "__version__", None) or getattr(mod, "VERSION", None)
    detail = f"{modname} {version}" if version else modname
    return ValidationIssue("OK", label, detail)


def _check_executable(name: str, label: str, *, required: bool) -> ValidationIssue:
    path = shutil.which(name)
    if path:
        return ValidationIssue("OK", label, path)
    level = "ERROR" if required else "WARN"
    return ValidationIssue(level, label, "not found on PATH")


def _check_optional_executable(name: str, label: str) -> ValidationIssue:
    path = shutil.which(name)
    if path:
        return ValidationIssue("OK", label, path)
    return ValidationIssue("OPTIONAL", label, "missing")


def _check_datapath(*, required: bool) -> ValidationIssue:
    val = os.environ.get("DATAPATH")
    if not val:
        level = "ERROR" if required else "WARN"
        return ValidationIssue(level, "RNAstructure DATAPATH", "not set")
    p = Path(val)
    if not p.exists() or not p.is_dir():
        level = "ERROR" if required else "WARN"
        return ValidationIssue(level, "RNAstructure DATAPATH", f"not a directory: {val}")
    return ValidationIssue("OK", "RNAstructure DATAPATH", val)


def _parse_args(args: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="mtrnafeat doctor")
    p.add_argument("--json", dest="json_path", type=Path, default=None,
                   help="Write machine-readable JSON to this path")
    return p.parse_args(args or [])


def _gather(cfg: Config, *, config_path: str | None) -> list[ValidationIssue]:
    rnastructure_required = cfg.fold_engine == "rnastructure"
    config_detail = config_path if config_path else "(default — no --config supplied)"
    issues = [
        _check_python(),
        _check_self(),
        _check_module("RNA", "ViennaRNA import"),
        _check_executable("RNAplfold", "RNAplfold availability", required=False),
        _check_executable("Fold", "RNAstructure Fold", required=rnastructure_required),
        _check_datapath(required=rnastructure_required),
        _check_optional_executable("DrTransformer", "DrTransformer (kinetic)"),
        ValidationIssue("OK", "Config load", config_detail),
        issues_writable(cfg.outdir),
        ValidationIssue("OK", "Working directory", os.getcwd()),
    ]
    return issues


def run(cfg: Config, args: list[str] | None = None) -> int:
    ns = _parse_args(args)
    config_path = os.environ.get("MTRNAFEAT_CONFIG_PATH")
    issues = _gather(cfg, config_path=config_path)

    print("mtrnafeat doctor")
    print("=" * 16)
    print(format_table(issues_to_rows(issues)))

    if ns.json_path:
        ns.json_path.parent.mkdir(parents=True, exist_ok=True)
        ns.json_path.write_text(json.dumps(serialize_issues(issues), indent=2))

    return 1 if has_error(issues) else 0
