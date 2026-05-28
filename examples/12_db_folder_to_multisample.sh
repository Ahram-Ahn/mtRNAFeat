#!/usr/bin/env bash
# Build a multi-sample config from every .db file in one folder, then run
# independent sample analysis. This is useful when the only input you have
# is a folder of DMS-derived .db files.
#
# Example:
#   DB_DIR=data/my_samples OUTDIR=runs/my_samples examples/12_db_folder_to_multisample.sh
#
# Set RUN=0 to only write the generated config.

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

DB_DIR="${DB_DIR:-data/my_samples}"
OUTDIR="${OUTDIR:-runs/$(date +%Y-%m-%d)-db-folder}"
CONFIG="${CONFIG:-$OUTDIR/generated.multisample.yaml}"
RUN="${RUN:-1}"

mkdir -p "$OUTDIR"

mapfile -t DBS < <(find "$DB_DIR" -maxdepth 1 -type f -name "*.db" | sort)
if [[ "${#DBS[@]}" -eq 0 ]]; then
  echo "[db-folder] no .db files found in $DB_DIR" >&2
  exit 1
fi

{
  echo "data_dir: $DB_DIR"
  echo "outdir: $OUTDIR"
  echo "db_files:"
  for db in "${DBS[@]}"; do
    base="$(basename "$db")"
    label="${base%.db}"
    label="$(printf '%s' "$label" | tr -c '[:alnum:]_' '_')"
    echo "  $label: $base"
  done
  echo "sample_annotation_species:"
  for db in "${DBS[@]}"; do
    base="$(basename "$db")"
    label="${base%.db}"
    label="$(printf '%s' "$label" | tr -c '[:alnum:]_' '_')"
    lower="$(printf '%s' "$label" | tr '[:upper:]' '[:lower:]')"
    if [[ "$lower" == *human* ]]; then
      echo "  $label: Human"
    elif [[ "$lower" == *yeast* ]]; then
      echo "  $label: Yeast"
    fi
  done
  cat <<'YAML'
seed: 42
sim_num_sequences: 500
gradient_steps: 11
gradient_seqs_per_step: 80
n_workers: 4
target_genes:
  - COX1
  - COX2
  - COX3
  - ATP86
  - CYTB
  - ND1
  - ND2
  - ND3
  - ND4L4
  - ND5
  - ND6
  - ATP9
  - COB
  - VAR1
fold_engine: vienna
YAML
} > "$CONFIG"

echo "[db-folder] wrote $CONFIG"
echo "[db-folder] samples:"
printf '  %s\n' "${DBS[@]}"

if [[ "$RUN" == "0" ]]; then
  echo "[db-folder] RUN=0, stopping before analysis."
  exit 0
fi

echo
echo "[db-folder] validating generated config"
"${MTRNAFEAT[@]}" validate-inputs --config "$CONFIG"

echo
echo "[db-folder] running independent sample analysis → $OUTDIR"
"${MTRNAFEAT[@]}" run-all --config "$CONFIG" --outdir "$OUTDIR" -- --parallel
