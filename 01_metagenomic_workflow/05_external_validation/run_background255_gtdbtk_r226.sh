#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep}
DONG=${DONG:-$ROOT/public_validation/Dong_NewbornSeep_PRJNA1313825_20260803}
PILOT=${PILOT:-$DONG/pilot4_v2}
WORK=${WORK:-$PILOT/11_background_MAG_enrichment_20260815}
MANIFEST=${MANIFEST:-$WORK/all_MQHQ_MAG_background_manifest.tsv}
GTDBTK=${GTDBTK:-/home/ps/anaconda3/envs/gtdbtk-2.4.1/bin/gtdbtk}
GTDB_ENV_BIN=${GTDB_ENV_BIN:-/home/ps/anaconda3/envs/gtdbtk-2.4.1/bin}
GTDB_RELEASE=${GTDB_RELEASE:-/home/ps/Research/databases/other_database/gtdbtk_226/release226}
CPUS=${CPUS:-48}
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
INPUTS=$WORK/00_gtdb_inputs
OUT=$WORK/01_gtdbtk_r226_all_MQHQ_retry1
LOGS=$WORK/logs

mkdir -p "$INPUTS" "$LOGS"
if [[ -s "$OUT/gtdbtk.bac120.summary.tsv" || -s "$OUT/gtdbtk.ar53.summary.tsv" ]]; then
  echo "GTDB-Tk summary already exists; refusing to overwrite: $OUT" >&2
  exit 2
fi
if [[ -d "$OUT" && -n "$(find "$OUT" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "Non-empty incomplete output exists; inspect before rerunning: $OUT" >&2
  exit 3
fi

python3 "$HERE/scripts/prepare_background_gtdb_inputs.py" \
  --manifest "$MANIFEST" \
  --outdir "$INPUTS" \
  > "$LOGS/background255_input_audit.log"

export GTDBTK_DATA_PATH="$GTDB_RELEASE"
export PATH="$GTDB_ENV_BIN:$PATH"
{
  date '+started=%F %T %z'
  "$GTDBTK" --version
  grep '^VERSION_DATA=' "$GTDB_RELEASE/metadata/metadata.txt"
  echo "GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"
  echo "command=$GTDBTK classify_wf --batchfile $INPUTS/all_MQHQ_MAG_GTDB_batchfile.tsv --out_dir $OUT --cpus $CPUS --pplacer_cpus $CPUS --skip_ani_screen --force"
} > "$LOGS/background255_gtdbtk_run_metadata.txt" 2>&1

"$GTDBTK" classify_wf \
  --batchfile "$INPUTS/all_MQHQ_MAG_GTDB_batchfile.tsv" \
  --out_dir "$OUT" \
  --cpus "$CPUS" \
  --pplacer_cpus "$CPUS" \
  --skip_ani_screen \
  --force

date '+completed=%F %T %z' >> "$LOGS/background255_gtdbtk_run_metadata.txt"
