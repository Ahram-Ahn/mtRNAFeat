"""`mtrnafeat validate-inputs` — pre-flight config and data validation.

Runs the full set of checks in `mtrnafeat.validation` against the user's
config and the data files it references. Prints a human-readable summary
table; with `--json PATH`, also writes a machine-readable list of
validation issues.

Exit status:
  - 0 if no ERRORs (warnings tolerated).
  - 1 if any required check fails.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from mtrnafeat.config import Config
from mtrnafeat.validation import (
    format_table,
    has_error,
    issues_to_rows,
    run_all_checks,
    serialize_issues,
)


def _parse_args(args: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="mtrnafeat validate-inputs")
    p.add_argument("--json", dest="json_path", type=Path, default=None,
                   help="Write machine-readable JSON to this path")
    return p.parse_args(args or [])


def run(cfg: Config, args: list[str] | None = None) -> int:
    ns = _parse_args(args)
    config_path = os.environ.get("MTRNAFEAT_CONFIG_PATH")
    if not config_path:
        print(
            "validate-inputs requires --config; pass it before the subcommand, e.g.\n"
            "    mtrnafeat validate-inputs --config configs/all.yaml",
            file=sys.stderr,
        )
        return 1

    issues = run_all_checks(config_path)

    print("Input validation summary")
    print("=" * 24)
    print(format_table(issues_to_rows(issues)))

    if ns.json_path:
        ns.json_path.parent.mkdir(parents=True, exist_ok=True)
        ns.json_path.write_text(json.dumps(serialize_issues(issues), indent=2))

    return 1 if has_error(issues) else 0
