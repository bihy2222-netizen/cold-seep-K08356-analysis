#!/usr/bin/env bash
set -euo pipefail

if [[ "${CONFIRM_JL_DNA_DOWNLOAD:-NO}" != "YES" ]]; then
  echo "Set CONFIRM_JL_DNA_DOWNLOAD=YES to download the paired JL_0.1 metagenome." >&2
  exit 2
fi

ROOT=/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra
RAW="$ROOT/01_reads/PRJNA831433_Jiaolong_DNA/JL_0.1_SRR19020591"
OUT="$ROOT/maggie_expression_validation_20260812/paired_DNA_screening/SRR19020591_JL_0.1_DNA"
LOG="$ROOT/logs/SRR19020591_JL_0.1_DNA"
HMM_DIR=/home/ps/ps1/data/bihongyu/cold_seep/Neighborhood_Analyses/hmm
THREADS=${THREADS:-16}
export PATH="/home/ps/anaconda3/envs/megahit/bin:/home/ps/anaconda3/bin:$PATH"
mkdir -p "$RAW" "$OUT" "$LOG"

for tool in curl md5sum gzip fastp megahit prodigal hmmsearch; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 3; }
done

available=$(df -PB1 "$ROOT" | awk 'NR==2 {print $4}')
required=$((3811089081 + 4110594514 + 40 * 1024 * 1024 * 1024))
if (( available < required )); then
  echo "Insufficient free space: available=$available required=$required" >&2
  exit 4
fi

download_one() {
  local label="$1" url="$2" expected_size="$3" expected_md5="$4" final="$5"
  local tmp="${final}.fresh.tmp"
  if [[ -e "$final" ]]; then
    local size md5
    size=$(stat -c %s "$final")
    md5=$(md5sum "$final" | awk '{print $1}')
    if [[ -s "$final" && "$size" == "$expected_size" && "$md5" == "$expected_md5" ]] && gzip -t "$final"; then
      printf '%s\texisting\tPASS\t%s\t%s\n' "$label" "$size" "$md5" >> "$LOG/download_validation.tsv"
      return 0
    fi
    echo "Existing final file failed validation and was not overwritten: $final" >&2
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
  echo "Failed after five fresh attempts: $label" >&2
  return 1
}

printf 'read_end\tsource\tvalidation\tactual_bytes\tactual_md5\n' > "$LOG/download_validation.tsv"
download_one R1 \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR190/091/SRR19020591/SRR19020591_1.fastq.gz \
  3811089081 bf934f629bca29cea51ead709123dc35 "$RAW/SRR19020591_1.fastq.gz"
download_one R2 \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR190/091/SRR19020591/SRR19020591_2.fastq.gz \
  4110594514 5ae44ea383f74179a3a2bb6578bf55fe "$RAW/SRR19020591_2.fastq.gz"

mkdir -p "$OUT/01_fastp" "$OUT/03_prodigal" "$OUT/04_hmm"
fastp \
  -i "$RAW/SRR19020591_1.fastq.gz" -I "$RAW/SRR19020591_2.fastq.gz" \
  -o "$OUT/01_fastp/SRR19020591.clean.R1.fastq.gz" \
  -O "$OUT/01_fastp/SRR19020591.clean.R2.fastq.gz" \
  --html "$OUT/01_fastp/SRR19020591.fastp.html" \
  --json "$OUT/01_fastp/SRR19020591.fastp.json" --thread "$THREADS"

megahit \
  -1 "$OUT/01_fastp/SRR19020591.clean.R1.fastq.gz" \
  -2 "$OUT/01_fastp/SRR19020591.clean.R2.fastq.gz" \
  --presets meta-sensitive --min-contig-len 300 -t "$THREADS" \
  -o "$OUT/02_megahit"

prodigal -p meta \
  -i "$OUT/02_megahit/final.contigs.fa" \
  -a "$OUT/03_prodigal/SRR19020591.proteins.faa" \
  -d "$OUT/03_prodigal/SRR19020591.CDS.fna" \
  -f gff -o "$OUT/03_prodigal/SRR19020591.genes.gff"

FAA="$OUT/03_prodigal/SRR19020591.proteins.faa"
hmmsearch --cpu "$THREADS" --tblout "$OUT/04_hmm/combined.all.tbl" \
  --domtblout "$OUT/04_hmm/combined.all.domtbl" \
  "$HMM_DIR/combined_iriA_aioA.hmm" "$FAA" > "$OUT/04_hmm/combined.all.txt"
hmmsearch --cpu "$THREADS" -T 640 --domT 640 \
  --tblout "$OUT/04_hmm/combined.T640.tbl" \
  --domtblout "$OUT/04_hmm/combined.T640.domtbl" \
  "$HMM_DIR/combined_iriA_aioA.hmm" "$FAA" > "$OUT/04_hmm/combined.T640.txt"
hmmsearch --cpu "$THREADS" --tblout "$OUT/04_hmm/iriA_new.all.tbl" \
  --domtblout "$OUT/04_hmm/iriA_new.all.domtbl" \
  "$HMM_DIR/iriA_new.hmm" "$FAA" > "$OUT/04_hmm/iriA_new.all.txt"
hmmsearch --cpu "$THREADS" --tblout "$OUT/04_hmm/aioA.all.tbl" \
  --domtblout "$OUT/04_hmm/aioA.all.domtbl" \
  "$HMM_DIR/aioA.hmm" "$FAA" > "$OUT/04_hmm/aioA.all.txt"

python3 "$(dirname "$0")/prepare_JL_DNA_candidate_review.py" \
  --combined "$OUT/04_hmm/combined.all.tbl" \
  --combined-strict "$OUT/04_hmm/combined.T640.tbl" \
  --iria "$OUT/04_hmm/iriA_new.all.tbl" \
  --aioa "$OUT/04_hmm/aioA.all.tbl" \
  --proteins "$OUT/03_prodigal/SRR19020591.proteins.faa" \
  --cds "$OUT/03_prodigal/SRR19020591.CDS.fna" \
  --gff "$OUT/03_prodigal/SRR19020591.genes.gff" \
  --outdir "$OUT/05_candidate_review"

cat > "$OUT/NEXT_STEP.txt" <<'EOF'
T640 denotes a high-scoring K08356/DMSOR-family candidate, not a functional IdrA call.
Review candidate proteins in the broad DMSOR tree and annotate their gene neighborhoods.
Only rows with final_class=strict_IdrA_associated, neighborhood_support=complete_DIRM_like,
and accept_for_personalized_reference=yes may be passed to build_personalized_idra_reference.py.
EOF
