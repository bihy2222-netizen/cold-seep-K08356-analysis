#!/usr/bin/env bash
set -euo pipefail

FAA_DIR="/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/all_50_10_dRep/dereplicated_genomes/bin_faa"
MODEL_DIR="/home/ps/ps1/data/bihongyu/cold_seep/Neighborhood_Analyses/hmm"
WORK_DIR="/home/ps/ps1/data/bihongyu/cold_seep/HMM_supplement/808MAG_IdrA_screen_20260807"
HMMSEARCH="/home/ps/mambaforge/bin/hmmsearch"
CPU="${CPU:-32}"

INPUT_DIR="$WORK_DIR/00_input"
MODEL_COPY_DIR="$WORK_DIR/01_models"
RESULT_DIR="$WORK_DIR/02_hmm_results"
TABLE_DIR="$WORK_DIR/03_tables"
LOG_DIR="$WORK_DIR/04_logs"
COMBINED_FAA="$INPUT_DIR/all_808MAG_proteins_prefixed.faa"

mkdir -p "$INPUT_DIR" "$MODEL_COPY_DIR" "$RESULT_DIR" "$TABLE_DIR" "$LOG_DIR"

mapfile -t faa_files < <(find "$FAA_DIR" -maxdepth 1 -type f -name '*.cds.faa' | sort)
if [[ "${#faa_files[@]}" -ne 808 ]]; then
  printf 'Expected 808 FAA files, found %s\n' "${#faa_files[@]}" >&2
  exit 1
fi

if find "$FAA_DIR" -maxdepth 1 -type f -name '*.cds.faa' -size 0 | grep -q .; then
  printf 'At least one FAA file is empty.\n' >&2
  exit 1
fi

cp "$MODEL_DIR/combined_iriA_aioA.hmm" "$MODEL_COPY_DIR/"
cp "$MODEL_DIR/iriA_new.hmm" "$MODEL_COPY_DIR/"
cp "$MODEL_DIR/aioA.hmm" "$MODEL_COPY_DIR/"

sha256sum "$MODEL_COPY_DIR"/*.hmm > "$TABLE_DIR/model_sha256.tsv"
printf 'model\trole\tthreshold\tinterpretation\n' > "$TABLE_DIR/model_roles.tsv"
printf 'combined_iriA_aioA.hmm\tprimary discovery screen\t640\tauthor repository README explicitly associates this model with threshold 640\n' >> "$TABLE_DIR/model_roles.tsv"
printf 'iriA_new.hmm\tIdrA-specific sensitivity screen\t640\texploratory because a model-specific 640 cutoff is not documented in the repository README\n' >> "$TABLE_DIR/model_roles.tsv"
printf 'aioA.hmm\tAioA comparator screen\t640\texploratory comparator\n' >> "$TABLE_DIR/model_roles.tsv"

if [[ ! -s "$COMBINED_FAA" ]]; then
  : > "$COMBINED_FAA"
  for faa in "${faa_files[@]}"; do
    mag="$(basename "$faa" .cds.faa)"
    awk -v mag="$mag" '/^>/ {sub(/^>/, ">" mag "|")} {print}' "$faa" >> "$COMBINED_FAA"
  done
fi

protein_count="$(grep -c '^>' "$COMBINED_FAA")"
duplicate_count="$(awk '/^>/{print $1}' "$COMBINED_FAA" | sort | uniq -d | wc -l)"
printf 'metric\tvalue\nFAA_files\t%s\nproteins\t%s\nduplicate_prefixed_IDs\t%s\n' \
  "${#faa_files[@]}" "$protein_count" "$duplicate_count" > "$TABLE_DIR/input_QC.tsv"

if [[ "$protein_count" -ne 1893843 || "$duplicate_count" -ne 0 ]]; then
  printf 'Combined protein QC failed: proteins=%s duplicates=%s\n' "$protein_count" "$duplicate_count" >&2
  exit 1
fi

run_search() {
  local label="$1"
  local model="$2"
  "$HMMSEARCH" --cpu "$CPU" -T 640 \
    --tblout "$RESULT_DIR/${label}_T640.tbl" \
    "$model" "$COMBINED_FAA" \
    > "$RESULT_DIR/${label}_T640.out"
}

run_search "combined_AioA_IdrA" "$MODEL_COPY_DIR/combined_iriA_aioA.hmm"
run_search "IdrA_specific" "$MODEL_COPY_DIR/iriA_new.hmm"
run_search "AioA_specific" "$MODEL_COPY_DIR/aioA.hmm"

for tbl in "$RESULT_DIR"/*_T640.tbl; do
  label="$(basename "$tbl" .tbl)"
  awk -v OFS='\t' -v model="$label" '!/^#/ {split($1,a,"|"); print model,a[1],a[2],$1,$5,$6}' "$tbl"
done | {
  printf 'screen\tMAG_ID\tprotein_ID\tprefixed_target_ID\tfull_sequence_Evalue\tfull_sequence_bitscore\n'
  cat
} > "$TABLE_DIR/all_T640_hits_long.tsv"

awk -F '\t' 'NR>1 {print $2}' "$TABLE_DIR/all_T640_hits_long.tsv" | sort -u > "$TABLE_DIR/T640_hit_MAG_IDs.txt"
awk -F '\t' 'NR>1 {print $4}' "$TABLE_DIR/all_T640_hits_long.tsv" | sort -u > "$TABLE_DIR/T640_hit_prefixed_protein_IDs.txt"

awk -F '\t' 'BEGIN{OFS="\t"} NR>1 {hits[$1]++; mags[$1 SUBSEP $2]=1; proteins[$1 SUBSEP $4]=1} END {print "screen","hits","unique_MAGs","unique_proteins"; for (s in hits) {m=0;p=0; for (k in mags) {split(k,a,SUBSEP); if(a[1]==s)m++} for(k in proteins){split(k,a,SUBSEP);if(a[1]==s)p++} print s,hits[s],m,p}}' \
  "$TABLE_DIR/all_T640_hits_long.tsv" | sort > "$TABLE_DIR/T640_screen_summary.tsv"

printf 'Completed %s\n' "$(date -Is)" > "$LOG_DIR/completed.txt"
