#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra
MAG="$ROOT/maggie_expression_validation_20260812"
SAMPLE=SRR19238834_JL_0.1
OUT="$MAG/maggie_diamond_blastx_pilot/$SAMPLE"
FAA="$MAG/competitive_DMSOR_plus_treeRefs_proteins.faa"
DB="$OUT/refs/competitive_DMSOR_plus_treeRefs_proteins"
R1="$ROOT/03_fastp/JL_0.1/SRR19238834.clean.R1.fastq.gz"
R2="$ROOT/03_fastp/JL_0.1/SRR19238834.clean.R2.fastq.gz"
mkdir -p "$OUT/refs" "$OUT/hits" "$OUT/logs" "$OUT/tables"

diamond makedb --in "$FAA" --db "$DB" > "$OUT/logs/diamond_makedb.log" 2>&1

for mate in R1 R2; do
  reads="$R1"; [[ "$mate" == R2 ]] && reads="$R2"
  diamond blastx \
    --db "$DB" --query "$reads" \
    --out "$OUT/hits/${SAMPLE}.${mate}.competitive63.blastx.tsv.gz" \
    --outfmt 6 qseqid qlen sseqid slen pident length mismatch gaps qstart qend sstart send evalue bitscore qframe \
    --very-sensitive --evalue 1e-5 --query-cover 60 --top 5 --threads 16 --compress 1 \
    > "$OUT/logs/${mate}.diamond.stdout.log" 2> "$OUT/logs/${mate}.diamond.stderr.log"
done

python3 - "$ROOT/03_fastp/JL_0.1/SRR19238834.fastp.json" "$OUT/total_reads_searched.tsv" <<'PY'
import json,sys
with open(sys.argv[1]) as handle:
    total=json.load(handle)["summary"]["after_filtering"]["total_reads"]
with open(sys.argv[2],"w") as handle:
    handle.write("sample\tread_end\ttotal_reads_searched\n")
    handle.write(f"SRR19238834_JL_0.1\tR1\t{total // 2}\n")
    handle.write(f"SRR19238834_JL_0.1\tR2\t{total // 2}\n")
PY

python3 "$MAG/scripts/summarize_maggie_diamond_blastx.py" \
  --hits-dir "$OUT/hits" --out-dir "$OUT/tables" --total-reads "$OUT/total_reads_searched.tsv"
