#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep}
DONG=${DONG:-$ROOT/public_validation/Dong_NewbornSeep_PRJNA1313825_20260803}
PILOT=${PILOT:-$DONG/pilot4_v2}
CLASSIFICATION=${CLASSIFICATION:-$PILOT/10_external82_terminal_classification_20260814_v3}
CONTEXT=${CONTEXT:-$PILOT/22_N30_16/10_maggie_mag_context_validation}
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

python3 "$HERE/scripts/summarize_dna_presence_control.py" \
  --classification "$CLASSIFICATION/external82_integrated_candidate_classification.tsv" \
  --neighborhoods "$CLASSIFICATION/external82_neighborhood_reviewed.tsv" \
  --bam "$PILOT/22_N30_16/05_read_mapping/22_N30_16.sorted.bam" \
  --sample 22_N30_16 \
  --outdir "$CONTEXT/06_per_base_depth/whole_assembly_primary_20260814"
