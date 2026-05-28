#!/usr/bin/env bash
# Same smoke run as 01, but with --parallel so every analysis stage fires
# concurrently as a subprocess. Useful for verifying that the parallel
# pipeline mode behaves correctly across all stages.

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

OUTDIR="${OUTDIR:-runs/smoke_parallel}"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "[smoke-parallel] parallel run → $OUTDIR"
"${MTRNAFEAT[@]}" run-all --config test_data/mini.config.yaml --outdir "$OUTDIR" \
       -- --parallel

echo
echo "[smoke-parallel] outputs:"
find "$OUTDIR" -maxdepth 2 -type f | sort
