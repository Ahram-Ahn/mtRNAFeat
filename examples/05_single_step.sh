#!/usr/bin/env bash
# Run a single subcommand. Useful for iterating on one analysis or one plot.
# Edit STEP below or pass it on the command line: STEP=nascent ./05_single_step.sh

set -euo pipefail
cd "$(dirname "$0")/.."

source .venv/bin/activate

STEP="${STEP:-landscape}"
CONFIG="${CONFIG:-test_data/mini.config.yaml}"
OUTDIR="${OUTDIR:-runs/single_${STEP}}"
EXTRA=("${@:1}")

mkdir -p "$OUTDIR"
echo "[single] $STEP → $OUTDIR (config=$CONFIG)"
mtrnafeat "$STEP" --config "$CONFIG" --outdir "$OUTDIR" "${EXTRA[@]}"

echo
echo "[single] outputs:"
find "$OUTDIR" -maxdepth 2 -type f | sort
