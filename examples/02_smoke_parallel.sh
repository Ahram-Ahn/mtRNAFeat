#!/usr/bin/env bash
# Same smoke run as 01, but with --parallel so every analysis stage fires
# concurrently as a subprocess. Useful for verifying that the parallel
# pipeline mode behaves correctly across all stages.

set -euo pipefail
cd "$(dirname "$0")/.."

source .venv/bin/activate

OUTDIR="${OUTDIR:-runs/smoke_parallel}"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "[smoke-parallel] parallel run → $OUTDIR"
mtrnafeat run-all --config test_data/mini.config.yaml --outdir "$OUTDIR" \
       -- --parallel --skip kinetic

echo
echo "[smoke-parallel] outputs:"
find "$OUTDIR" -maxdepth 2 -type f | sort
