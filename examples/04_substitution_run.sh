#!/usr/bin/env bash
# Run only the substitution-thermo permutation test.
# For each (species, gene), the wild-type CDS is compared to codon-aware
# null pools under plain Vienna MFE. CoFold is handled by the separate
# `mtrnafeat cofold` stage.

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

CONFIG="${CONFIG:-test_data/mini.config.yaml}"
OUTDIR="${OUTDIR:-runs/substitution}"
N="${N:-200}"
MAX_NT="${MAX_NT:-300}"

mkdir -p "$OUTDIR"
echo "[substitution] running with N=$N max_nt=$MAX_NT → $OUTDIR"

"${MTRNAFEAT[@]}" substitution --config "$CONFIG" --outdir "$OUTDIR" -- --n "$N" --max-nt "$MAX_NT"

echo
echo "[substitution] outputs:"
find "$OUTDIR" -maxdepth 2 -type f | sort
echo
echo "[substitution] summary table:"
column -t -s, < "$OUTDIR/tables/substitution_thermo_summary.csv" | head -25
