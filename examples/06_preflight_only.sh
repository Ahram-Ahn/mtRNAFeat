#!/usr/bin/env bash
# Check the local environment and input files before spending time on a
# full analysis. This does not create plots or expensive simulations.

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

CONFIG="${CONFIG:-configs/all.yaml}"

echo "[preflight] doctor: $CONFIG"
"${MTRNAFEAT[@]}" doctor --config "$CONFIG"

echo
echo "[preflight] validate inputs: $CONFIG"
"${MTRNAFEAT[@]}" validate-inputs --config "$CONFIG"
