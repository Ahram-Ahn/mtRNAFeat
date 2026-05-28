#!/usr/bin/env bash
# Re-run only landscape/ and features/. Use this when you are checking GC
# nulls, experimental overlays, structural motifs, phase-space contours, or
# span ECDFs without rerunning every gene-level stage.

set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}src"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/mtrnafeat-mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/mtrnafeat-cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
MTRNAFEAT=("$PYTHON" -m mtrnafeat.cli)

CONFIG="${CONFIG:-configs/all.yaml}"
OUTDIR="${OUTDIR:-runs/$(date +%Y-%m-%d)-landscape-features}"

mkdir -p "$OUTDIR"

echo "[landscape-features] landscape → $OUTDIR"
"${MTRNAFEAT[@]}" landscape --config "$CONFIG" --outdir "$OUTDIR"

echo
echo "[landscape-features] features → $OUTDIR"
"${MTRNAFEAT[@]}" features --config "$CONFIG" --outdir "$OUTDIR"

echo
echo "[landscape-features] outputs:"
find "$OUTDIR" -maxdepth 2 \( -path "*/landscape/*" -o -path "*/features/*" -o -path "*/tables/*" \) -type f | sort
