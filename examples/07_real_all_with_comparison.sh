#!/usr/bin/env bash
# Full built-in Human + Yeast run with comparison stages explicitly kept.
# Use this for the canonical two-sample analysis where compare/ and
# substitution/ are meaningful.

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

CONFIG="${CONFIG:-configs/all.yaml}"
OUTDIR="${OUTDIR:-runs/$(date +%Y-%m-%d)-all-with-comparison}"

mkdir -p "$OUTDIR"

echo "[real-all-comparison] $CONFIG → $OUTDIR"
"${MTRNAFEAT[@]}" run-all --config "$CONFIG" --outdir "$OUTDIR" -- \
  --parallel \
  --include-comparison

echo
echo "[real-all-comparison] top-level outputs:"
find "$OUTDIR" -maxdepth 2 -type f | sort
