#!/usr/bin/env bash
# Run a single subcommand. Useful for iterating on one analysis or one plot.
# Edit STEP below or pass it on the command line:
#   STEP=tis ./05_single_step.sh
#   STEP=local-probability CONFIG=configs/all.yaml ./05_single_step.sh
# Pass subcommand-specific args after a literal `--`:
#   STEP=substitution ./05_single_step.sh -- --n 50 --max-nt 300
#
# Valid STEP values match `mtrnafeat <subcommand>`:
#   stats, landscape, features, window, local-probability, significance,
#   tis, compare, substitution, cofold, gene-panel, kinetic, plot.

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
