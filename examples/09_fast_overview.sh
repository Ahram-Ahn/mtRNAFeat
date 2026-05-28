#!/usr/bin/env bash
# Quick overview run for a new config. It keeps the fast summary stages and
# skips the heavier per-gene/probability/cofold stages.
#
# Runs: stats, landscape, features, tis, gene_panel
# Skips: window, local_probability, structure_deviation, cofold, compare,
# substitution

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

CONFIG="${CONFIG:-configs/all.yaml}"
OUTDIR="${OUTDIR:-runs/$(date +%Y-%m-%d)-fast-overview}"

mkdir -p "$OUTDIR"

echo "[fast-overview] $CONFIG → $OUTDIR"
"${MTRNAFEAT[@]}" run-all --config "$CONFIG" --outdir "$OUTDIR" -- \
  --parallel \
  --skip window,local_probability,structure_deviation,cofold,compare,substitution

echo
echo "[fast-overview] key outputs:"
find "$OUTDIR" -maxdepth 2 \( -name "*.csv" -o -name "*.svg" \) | sort
