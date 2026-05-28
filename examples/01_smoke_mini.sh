#!/usr/bin/env bash
# Tiny smoke run on the test_data/ fixtures (ND6 + ATP9, ~1.5 kb total).
# Expected wall time: <2 minutes on a laptop. Useful for iterating on plot
# code or for verifying a fresh checkout works end-to-end.

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

OUTDIR="${OUTDIR:-runs/smoke}"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "[smoke] sequential run → $OUTDIR"
"${MTRNAFEAT[@]}" run-all --config test_data/mini.config.yaml --outdir "$OUTDIR"
echo
echo "[smoke] outputs:"
find "$OUTDIR" -maxdepth 2 -type f | sort
