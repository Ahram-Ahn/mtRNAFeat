#!/usr/bin/env bash
# Create a temporary one-gene config from an existing config and run the
# sample pipeline. This is useful for debugging plots or folding settings
# on COX1/ND6/etc. without processing every transcript.
#
# Examples:
#   GENES=COX1 examples/11_one_gene_debug.sh
#   GENES=ND6 CONFIG=test_data/mini.config.yaml examples/11_one_gene_debug.sh
#   GENES=COX1,COX2 SKIP=cofold examples/11_one_gene_debug.sh

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

CONFIG="${CONFIG:-configs/all.yaml}"
GENES="${GENES:-COX1}"
OUTDIR="${OUTDIR:-runs/debug_$(echo "$GENES" | tr ',' '_')}"
SKIP="${SKIP:-cofold,compare,substitution}"
TMP_CONFIG="$OUTDIR/debug.config.yaml"

mkdir -p "$OUTDIR"

python - "$CONFIG" "$TMP_CONFIG" "$GENES" "$OUTDIR" <<'PY'
from pathlib import Path
import sys
import yaml

src, dst, genes, outdir = sys.argv[1:5]
with open(src) as handle:
    cfg = yaml.safe_load(handle)
cfg["target_genes"] = [g.strip() for g in genes.split(",") if g.strip()]
cfg["outdir"] = outdir

# Keep debug runs quick. Increase these in the generated config if you need
# publication-quality simulation clouds.
cfg["sim_num_sequences"] = min(int(cfg.get("sim_num_sequences", 1500)), 200)
cfg["gradient_seqs_per_step"] = min(int(cfg.get("gradient_seqs_per_step", 200)), 50)

Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

ARGS=(--parallel)
if [[ -n "$SKIP" ]]; then
  ARGS+=(--skip "$SKIP")
fi

echo "[one-gene-debug] genes=$GENES"
echo "[one-gene-debug] generated config=$TMP_CONFIG"
echo "[one-gene-debug] outdir=$OUTDIR"
"${MTRNAFEAT[@]}" run-all --config "$TMP_CONFIG" --outdir "$OUTDIR" -- "${ARGS[@]}"

echo
echo "[one-gene-debug] outputs:"
find "$OUTDIR" -maxdepth 2 -type f | sort
