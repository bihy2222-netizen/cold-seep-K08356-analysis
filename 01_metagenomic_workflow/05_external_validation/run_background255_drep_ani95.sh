#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep}
DONG=${DONG:-$ROOT/public_validation/Dong_NewbornSeep_PRJNA1313825_20260803}
PILOT=${PILOT:-$DONG/pilot4_v2}
WORK=${WORK:-$PILOT/11_background_MAG_enrichment_20260815}
INPUTS=${INPUTS:-$WORK/00_gtdb_inputs}
DREP_ENV_BIN=${DREP_ENV_BIN:-/home/ps/anaconda3/envs/derep/bin}
DREP=${DREP:-$DREP_ENV_BIN/dRep}
CPUS=${CPUS:-32}
OUT=${OUT:-$WORK/03_drep_ANI95_AF30}
LOGS=$WORK/logs
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

mkdir -p "$INPUTS" "$LOGS"
python3 "$HERE/scripts/prepare_background_gtdb_inputs.py" \
  --manifest "$WORK/all_MQHQ_MAG_background_manifest.tsv" \
  --outdir "$INPUTS" \
  > "$LOGS/background255_input_audit.drep.log"

if [[ -d "$OUT" && -n "$(find "$OUT" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "Non-empty dRep output exists; refusing to overwrite: $OUT" >&2
  exit 2
fi

export PATH="$DREP_ENV_BIN:$PATH"
{
  date '+started=%F %T %z'
  echo "command=$DREP dereplicate $OUT -g $INPUTS/drep_genome_paths.txt --genomeInfo $INPUTS/drep_genomeInfo.csv -p $CPUS -comp 50 -con 10 -pa 0.90 -sa 0.95 -nc 0.30 --S_algorithm fastANI --skip_plots"
  echo "definition=ANI95 secondary clusters with minimum aligned fraction 0.30; representatives will be re-audited lexicographically"
} > "$LOGS/background255_drep_ANI95_run_metadata.txt"

"$DREP" dereplicate "$OUT" \
  -g "$INPUTS/drep_genome_paths.txt" \
  --genomeInfo "$INPUTS/drep_genomeInfo.csv" \
  -p "$CPUS" \
  -comp 50 \
  -con 10 \
  -pa 0.90 \
  -sa 0.95 \
  -nc 0.30 \
  --S_algorithm fastANI \
  --skip_plots

date '+completed=%F %T %z' >> "$LOGS/background255_drep_ANI95_run_metadata.txt"
