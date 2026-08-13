#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra
MAG="$ROOT/maggie_expression_validation_20260812"
DNA="$MAG/paired_DNA_screening/SRR19020591_JL_0.1_DNA"
RNA="$ROOT/03_fastp/JL_0.1"
PERSONAL="$DNA/06_personalized_reference"
SCRIPTS="$MAG/scripts"
SAMPLE=SRR19238834_JL_0.1_personalized

REVIEW="$DNA/05_candidate_review/JL_0.1_K08356_candidate_review.tsv"
CDS="$DNA/03_prodigal/SRR19020591.CDS.fna"
FAA="$DNA/03_prodigal/SRR19020591.proteins.faa"
CLEAN_R1="$RNA/SRR19238834.clean.R1.fastq.gz"
CLEAN_R2="$RNA/SRR19238834.clean.R2.fastq.gz"
COMP63="$MAG/competitive_DMSOR_plus_treeRefs_CDS.fna"

for path in "$REVIEW" "$CDS" "$FAA" "$CLEAN_R1" "$CLEAN_R2" "$COMP63"; do
  [[ -s "$path" ]] || { echo "Required input missing: $path" >&2; exit 2; }
done

python3 "$SCRIPTS/build_personalized_idra_reference.py" \
  --review-tsv "$REVIEW" --cds "$CDS" --proteins "$FAA" --outdir "$PERSONAL"

PRIMARY="$PERSONAL/JL_0.1_personalized_strict_IdrA_CDS.fna"
COMP_PERSONAL="$PERSONAL/JL_0.1_personalized_plus_competitive63_CDS.fna"
cat "$PRIMARY" "$COMP63" > "$COMP_PERSONAL"

python3 "$SCRIPTS/run_maggie_mapping_pilot.py" \
  --root "$MAG" --sample "$SAMPLE" --r1 "$CLEAN_R1" --r2 "$CLEAN_R2" \
  --primary-ref "$PRIMARY" --competitive-ref "$COMP_PERSONAL" \
  --primary-mode personalized_strict_IdrA \
  --competitive-mode personalized_plus_competitive63 --threads "${THREADS:-16}"

