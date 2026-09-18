#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep/public_validation/Dong_NewbornSeep_PRJNA1313825_20260803}"
WORK="$ROOT/pilot4_v2/11_background_MAG_enrichment_20260815"
SCRIPT_DIR="$ROOT/scripts/external82_terminal_classification"
BIN_STATS="$WORK/02_bin_level_enrichment"
DREP="$WORK/03_drep_ANI95_AF30"
OUT="$WORK/04_ANI95_enrichment"

mkdir -p "$OUT" "$WORK/logs"

while [[ ! -s "$BIN_STATS/01_all_MQHQ_MAG_GTDB_taxonomy.tsv" ]]; do
    sleep 60
done

python3 "$SCRIPT_DIR/calculate_ani95_host_enrichment.py" \
    --background-taxonomy "$BIN_STATS/01_all_MQHQ_MAG_GTDB_taxonomy.tsv" \
    --cdb "$DREP/data_tables/Cdb.csv" \
    --outdir "$OUT"

sha256sum "$OUT"/*.tsv > "$OUT/SHA256SUMS.txt"
date '+ANI95 enrichment complete: %F %T %Z'
