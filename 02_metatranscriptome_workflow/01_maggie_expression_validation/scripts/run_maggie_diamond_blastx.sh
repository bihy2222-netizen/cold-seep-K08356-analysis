#!/usr/bin/env bash
set -euo pipefail

MAG=/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra/maggie_expression_validation_20260812
OUT="$MAG/maggie_diamond_blastx_pilot"
DB="$OUT/refs/competitive_DMSOR_plus_treeRefs_proteins"
FAA="$MAG/competitive_DMSOR_plus_treeRefs_proteins.faa"
THREADS=16

mkdir -p "$OUT/refs" "$OUT/logs" "$OUT/hits" "$OUT/tables"

if [[ ! -s "${DB}.dmnd" ]]; then
  diamond makedb --in "$FAA" --db "$DB" \
    > "$OUT/logs/diamond_makedb.log" 2>&1
fi

run_one() {
  local sample="$1"
  local mate="$2"
  local reads="$3"
  local hits="$OUT/hits/${sample}.${mate}.competitive63.blastx.tsv.gz"

  [[ -s "$reads" ]] || { echo "Missing reads: $reads" >&2; return 1; }
  [[ -s "$hits" ]] && { echo "Keeping existing result: $hits"; return 0; }

  diamond blastx \
    --db "$DB" \
    --query "$reads" \
    --out "$hits" \
    --outfmt 6 qseqid qlen sseqid slen pident length mismatch gaps qstart qend sstart send evalue bitscore qframe \
    --very-sensitive \
    --evalue 1e-5 \
    --query-cover 60 \
    --top 5 \
    --threads "$THREADS" \
    --compress 1 \
    > "$OUT/logs/${sample}.${mate}.diamond.stdout.log" \
    2> "$OUT/logs/${sample}.${mate}.diamond.stderr.log"
}

run_one SRR16610255_QDN-W01B-1949 R1 \
  /home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra/03_fastp/QDN-W01B-1949/SRR16610255.clean.R1.fastq.gz
run_one SRR16610255_QDN-W01B-1949 R2 \
  /home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra/03_fastp/QDN-W01B-1949/SRR16610255.clean.R2.fastq.gz
run_one SRR16610253_QDN-W04B-4900 R1 \
  /home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra/03_fastp/QDN-W04B-4900/SRR16610253.clean.R1.fastq.gz
run_one SRR16610253_QDN-W04B-4900 R2 \
  /home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra/03_fastp/QDN-W04B-4900/SRR16610253.clean.R2.fastq.gz

python3 "$MAG/scripts/summarize_maggie_diamond_blastx.py" \
  --hits-dir "$OUT/hits" \
  --out-dir "$OUT/tables" \
  --total-reads "$OUT/total_reads_searched.tsv"
