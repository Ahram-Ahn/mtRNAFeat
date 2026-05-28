#!/usr/bin/env bash
# Full pipeline on both species (Human + Yeast .db files) in one run.
# Covers all 14 unique mt-mRNAs (4 shared, 7 human-only, 3 yeast-only).
# Expected wall time: 30–60 min in parallel mode.

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

OUTDIR="${OUTDIR:-runs/$(date +%Y-%m-%d)-all}"
mkdir -p "$OUTDIR"

echo "[real-all] parallel run → $OUTDIR"
"${MTRNAFEAT[@]}" run-all --config configs/all.yaml --outdir "$OUTDIR" -- --parallel
echo
echo "[real-all] top-level CSVs:"
find "$OUTDIR" -maxdepth 2 -name "*.csv" | sort
