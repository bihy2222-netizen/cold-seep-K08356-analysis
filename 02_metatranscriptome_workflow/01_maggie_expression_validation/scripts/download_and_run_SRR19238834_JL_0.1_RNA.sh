#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra
MAG="$ROOT/maggie_expression_validation_20260812"
RAW="$ROOT/01_reads/PRJNA831433_Jiaolong_RNA/JL_0.1_SRR19238834"
FASTP_OUT="$ROOT/03_fastp/JL_0.1"
SAMPLE=SRR19238834_JL_0.1
LOG="$ROOT/logs/SRR19238834_JL_0.1_RNA"
mkdir -p "$RAW" "$FASTP_OUT" "$LOG"

download_one() {
  local label="$1" url="$2" expected_size="$3" expected_md5="$4" final="$5"
  local tmp="${final}.fresh.tmp"

  if [[ -s "$final" ]]; then
    local size md5
    size=$(stat -c %s "$final")
    md5=$(md5sum "$final" | awk '{print $1}')
    if [[ "$size" == "$expected_size" && "$md5" == "$expected_md5" ]] && gzip -t "$final"; then
      printf '%s\texisting\tPASS\t%s\t%s\n' "$label" "$size" "$md5" >> "$LOG/download_validation.tsv"
      return 0
    fi
    printf 'Existing file failed validation and was not overwritten: %s\n' "$final" >&2
    return 1
  fi

  for attempt in 1 2 3 4 5; do
    rm -f "$tmp"
    printf '%s attempt %s started %s\n' "$label" "$attempt" "$(date -Is)" >> "$LOG/download_attempts.log"
    if curl --fail --location --retry 0 --output "$tmp" "$url" \
      >> "$LOG/${label}.curl.stdout.log" 2>> "$LOG/${label}.curl.stderr.log"; then
      local size md5
      size=$(stat -c %s "$tmp")
      md5=$(md5sum "$tmp" | awk '{print $1}')
      if [[ "$size" == "$expected_size" && "$md5" == "$expected_md5" ]] && gzip -t "$tmp"; then
        mv "$tmp" "$final"
        printf '%s\tdownloaded_attempt_%s\tPASS\t%s\t%s\n' "$label" "$attempt" "$size" "$md5" >> "$LOG/download_validation.tsv"
        return 0
      fi
      printf '%s attempt %s validation failed size=%s md5=%s\n' "$label" "$attempt" "$size" "$md5" >> "$LOG/download_attempts.log"
    fi
    rm -f "$tmp"
  done
  printf 'Failed after 5 fresh attempts: %s\n' "$label" >&2
  return 1
}

printf 'read_end\tsource\tvalidation\tactual_bytes\tactual_md5\n' > "$LOG/download_validation.tsv"
download_one R1 \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/034/SRR19238834/SRR19238834_1.fastq.gz \
  533274148 0fda5686a6b68aebb10fff290568798d "$RAW/SRR19238834_1.fastq.gz"
download_one R2 \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR192/034/SRR19238834/SRR19238834_2.fastq.gz \
  600860868 bca0b0545725ea287a441bdde2caf546 "$RAW/SRR19238834_2.fastq.gz"

fastp \
  -i "$RAW/SRR19238834_1.fastq.gz" \
  -I "$RAW/SRR19238834_2.fastq.gz" \
  -o "$FASTP_OUT/SRR19238834.clean.R1.fastq.gz" \
  -O "$FASTP_OUT/SRR19238834.clean.R2.fastq.gz" \
  --html "$FASTP_OUT/SRR19238834.fastp.html" \
  --json "$FASTP_OUT/SRR19238834.fastp.json" \
  --thread 16 \
  > "$FASTP_OUT/SRR19238834.fastp.stdout.log" \
  2> "$FASTP_OUT/SRR19238834.fastp.stderr.log"

python3 "$MAG/scripts/run_maggie_mapping_pilot.py" \
  --root "$MAG" \
  --sample "$SAMPLE" \
  --r1 "$FASTP_OUT/SRR19238834.clean.R1.fastq.gz" \
  --r2 "$FASTP_OUT/SRR19238834.clean.R2.fastq.gz" \
  --primary-ref "$MAG/high_confidence_IdrA_CDS.fna" \
  --competitive-ref "$MAG/competitive_DMSOR_plus_treeRefs_CDS.fna" \
  --threads 16

bash "$MAG/scripts/run_SRR19238834_JL_0.1_diamond_blastx.sh"
