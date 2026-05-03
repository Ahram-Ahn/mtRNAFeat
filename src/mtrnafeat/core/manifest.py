"""Per-run reproducibility manifest.

Each ``mtrnafeat`` invocation writes ``<outdir>/run_manifest.json`` with
the package version, the resolved upstream-tool versions (ViennaRNA,
RNAstructure when present), the command and its full argv, the active
git commit (if running in a git checkout), the seed, the config path,
and an ISO-8601 UTC timestamp.

Reviewers and downstream consumers can use this manifest to (a) reproduce
a run with the same tool stack and (b) audit drift if outputs change
between runs (e.g. a ViennaRNA upgrade silently shifts pair probabilities).

The manifest is overwritten on each invocation; for ``run-all`` in
parallel mode the last subprocess to finish wins, but all subprocesses
share the same version stack so this is safe.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from mtrnafeat import __version__ as _MTRNAFEAT_VERSION


def _viennarna_version() -> str | None:
    try:
        import RNA  # type: ignore
    except ImportError:
        return None
    for attr in ("__version__", "VERSION"):
        v = getattr(RNA, attr, None)
        if v:
            return str(v)
    return "unknown (loaded)"


def _rnastructure_version() -> str | None:
    """Best-effort RNAstructure version probe.

    RNAstructure binaries don't have a stable ``--version`` flag across
    releases; we read ``$DATAPATH`` and the ``Fold`` binary location as
    a coarse fingerprint, which is what reviewers need to identify the
    install.
    """
    fold_path = shutil.which("Fold")
    if fold_path is None:
        return None
    datapath = os.environ.get("DATAPATH", "")
    # Some recent RNAstructure builds support ``--version``; older ones
    # treat it as an illegal option and emit a usage message. Accept the
    # version string only when stdout looks like a single short version
    # line; otherwise fall back to the binary-path fingerprint, which is
    # what reviewers actually need to identify the install.
    try:
        proc = subprocess.run(
            [fold_path, "--version"],
            capture_output=True, text=True, timeout=5,
        )
        if proc.returncode == 0 and proc.stdout.strip():
            first_line = proc.stdout.strip().splitlines()[0]
            if len(first_line) <= 80 and "illegal" not in first_line.lower():
                return f"{first_line} (binary={fold_path}, DATAPATH={datapath})"
    except Exception:
        pass
    return f"binary={fold_path}, DATAPATH={datapath}"


def _git_commit(repo_root: Path) -> str | None:
    if shutil.which("git") is None:
        return None
    try:
        proc = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        if proc.returncode == 0:
            sha = proc.stdout.strip()
            # Mark dirty trees so the manifest can't be misread as
            # reflecting a clean commit when local edits are present.
            dirty = subprocess.run(
                ["git", "-C", str(repo_root), "status", "--porcelain"],
                capture_output=True, text=True, timeout=5,
            )
            if dirty.returncode == 0 and dirty.stdout.strip():
                return f"{sha} (dirty)"
            return sha
    except Exception:
        return None
    return None


def write_run_manifest(
    outdir: Path,
    *,
    command: str,
    argv: list[str],
    config_path: Path | None,
    seed: int | None,
) -> Path:
    """Serialize the run manifest to ``<outdir>/run_manifest.json``."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    repo_root = Path(__file__).resolve().parents[3]
    manifest: dict[str, Any] = {
        "mtrnafeat_version": _MTRNAFEAT_VERSION,
        "command": command,
        "argv": list(argv),
        "config_path": str(config_path) if config_path is not None else None,
        "seed": seed,
        "viennarna_version": _viennarna_version(),
        "rnastructure_version": _rnastructure_version(),
        "git_commit": _git_commit(repo_root),
        "python_version": sys.version.split()[0],
        "platform": sys.platform,
        "timestamp_utc": _dt.datetime.now(_dt.UTC).isoformat(
            timespec="seconds"
        ),
    }
    path = outdir / "run_manifest.json"
    with open(path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return path
