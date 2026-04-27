#!/usr/bin/env bash
# Run only the substitution-thermo permutation test.
# This is the most novel piece of the package — for each (species, gene),
# the wild-type CDS is compared to N synonymous-recoding null pools under CoFold.

set -euo pipefail
cd "$(dirname "$0")/.."

source .venv/bin/activate

CONFIG="${CONFIG:-test_data/mini.config.yaml}"
OUTDIR="${OUTDIR:-runs/substitution}"
N="${N:-200}"
MAX_NT="${MAX_NT:-300}"

mkdir -p "$OUTDIR"
echo "[substitution] running with N=$N max_nt=$MAX_NT → $OUTDIR"

mtrnafeat substitution --config "$CONFIG" --outdir "$OUTDIR" -- --n "$N" --max-nt "$MAX_NT"

echo
echo "[substitution] outputs:"
find "$OUTDIR" -maxdepth 2 -type f | sort
echo
echo "[substitution] summary table:"
column -t -s, < "$OUTDIR/tables/substitution_thermo_summary.csv" | head -25
