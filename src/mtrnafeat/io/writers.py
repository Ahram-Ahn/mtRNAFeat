"""Canonical CSV writer: deterministic columns, fixed float format, atomic write."""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pandas as pd

from mtrnafeat.constants import CSV_FLOAT_FORMAT


def tables_csv(df: pd.DataFrame, out_root, basename: str) -> Path:
    """Write `df` to `<out_root>/tables/<basename>.csv` with canonical formatting.

    Centralizes all CSV-only outputs (stats, tis, concordance summaries,
    substitution summaries, comparative tables) so they're easy to find as a
    flat collection of analysis results — instead of being scattered across
    one subdirectory per stage.
    """
    out_dir = Path(out_root) / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    return canonical_csv(df, out_dir / f"{basename}.csv")


def canonical_csv(df: pd.DataFrame, path: str | Path) -> Path:
    """Write `df` to `path` with stable column order and float formatting.

    - Columns sorted alphabetically except 'Gene', 'Species', 'Condition' (if
      present) come first.
    - Floats written with CSV_FLOAT_FORMAT.
    - Atomic write (tmp file in same dir, then os.replace).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    pinned = [c for c in ("Gene", "Species", "Condition", "Transcript", "Transcript_Length",
                          "Window_Index", "Window_Start_1based", "Window_End_1based") if c in df.columns]
    rest = sorted(c for c in df.columns if c not in pinned)
    out = df[pinned + rest]

    fd, tmp_name = tempfile.mkstemp(suffix=".csv", dir=str(path.parent))
    os.close(fd)
    try:
        out.to_csv(tmp_name, index=False, float_format=CSV_FLOAT_FORMAT)
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.remove(tmp_name)
        raise
    return path
