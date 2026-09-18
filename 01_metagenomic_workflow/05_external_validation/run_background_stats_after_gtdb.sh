#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep}
DONG=${DONG:-$ROOT/public_validation/Dong_NewbornSeep_PRJNA1313825_20260803}
PILOT=${PILOT:-$DONG/pilot4_v2}
WORK=${WORK:-$PILOT/11_background_MAG_enrichment_20260815}
GTDB_OUT=${GTDB_OUT:-$WORK/01_gtdbtk_r226_all_MQHQ_retry1}
INPUTS=${INPUTS:-$WORK/00_gtdb_inputs}
CLASSIFICATION=${CLASSIFICATION:-$PILOT/10_external82_terminal_classification_20260814_v3/external82_integrated_candidate_classification.tsv}
STATS_OUT=${STATS_OUT:-$WORK/02_bin_level_enrichment}
GTDB_PID=${GTDB_PID:-3719577}
MAX_POLLS=${MAX_POLLS:-144}
POLL_SECONDS=${POLL_SECONDS:-300}
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

BAC=$GTDB_OUT/classify/gtdbtk.bac120.summary.tsv
ARC=$GTDB_OUT/classify/gtdbtk.ar53.summary.tsv
for ((poll=1; poll<=MAX_POLLS; poll++)); do
  if ! kill -0 "$GTDB_PID" 2>/dev/null; then
    break
  fi
  sleep "$POLL_SECONDS"
done
if kill -0 "$GTDB_PID" 2>/dev/null; then
  echo "Timed out while GTDB-Tk PID was still running: $GTDB_PID" >&2
  exit 2
fi
if [[ ! -s "$BAC" && ! -s "$ARC" ]]; then
  echo "GTDB-Tk ended without classification summaries" >&2
  exit 3
fi

mkdir -p "$STATS_OUT"
python3 "$HERE/scripts/calculate_background_host_enrichment.py" \
  --input-crosswalk "$INPUTS/all_MQHQ_MAG_GTDB_input_crosswalk.tsv" \
  --gtdb-summary "$BAC" "$ARC" \
  --candidate-classification "$CLASSIFICATION" \
  --outdir "$STATS_OUT"

find "$INPUTS" "$GTDB_OUT" "$STATS_OUT" -type f -print0 \
  | sort -z \
  | xargs -0 sha256sum > "$WORK/12_checksums.sha256"
