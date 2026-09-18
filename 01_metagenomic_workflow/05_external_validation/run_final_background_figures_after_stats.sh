#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep/public_validation/Dong_NewbornSeep_PRJNA1313825_20260803}"
WORK="$ROOT/pilot4_v2/11_background_MAG_enrichment_20260815"
SCRIPT_DIR="$ROOT/scripts/external82_terminal_classification"
BIN_STATS="$WORK/02_bin_level_enrichment"
ANI_STATS="$WORK/04_ANI95_enrichment"
OUT="$WORK/05_figures"

while [[ ! -s "$BIN_STATS/04_strict_Fisher_OR_CI.tsv" || \
         ! -s "$ANI_STATS/07_strict_dereplicated_Fisher_OR_CI.tsv" ]]; do
    sleep 60
done

mkdir -p "$OUT"
python3 "$SCRIPT_DIR/plot_background_host_enrichment.py" \
    --bin-dir "$BIN_STATS" \
    --ani-dir "$ANI_STATS" \
    --outdir "$OUT"

sha256sum "$BIN_STATS"/* "$ANI_STATS"/* "$OUT"/* > "$WORK/12_checksums.sha256"
date '+Final background tables and figures complete: %F %T %Z'
