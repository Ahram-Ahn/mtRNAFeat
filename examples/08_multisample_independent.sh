#!/usr/bin/env bash
# Independent analysis of many samples. Use this when your config has more
# than the canonical Human/Yeast pair, for example four yeast KOs and two
# human treatment samples. run-all will skip compare/ and substitution/
# unless you force --include-comparison.

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

CONFIG="${CONFIG:-configs/my-experiment.yaml}"
OUTDIR="${OUTDIR:-runs/$(date +%Y-%m-%d)-multisample}"
SKIP="${SKIP:-}"

mkdir -p "$OUTDIR"

ARGS=(--parallel)
if [[ -n "$SKIP" ]]; then
  ARGS+=(--skip "$SKIP")
fi

echo "[multisample] independent sample run"
echo "[multisample] config=$CONFIG"
echo "[multisample] outdir=$OUTDIR"
if [[ -n "$SKIP" ]]; then
  echo "[multisample] skipping=$SKIP"
fi

"${MTRNAFEAT[@]}" run-all --config "$CONFIG" --outdir "$OUTDIR" -- "${ARGS[@]}"

echo
echo "[multisample] outputs:"
find "$OUTDIR" -maxdepth 2 -type f | sort
