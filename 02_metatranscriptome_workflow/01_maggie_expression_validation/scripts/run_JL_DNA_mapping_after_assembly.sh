#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra
DNA="$ROOT/maggie_expression_validation_20260812/paired_DNA_screening/SRR19020591_JL_0.1_DNA"
ASSEMBLY="$DNA/02_megahit/final.contigs.fa"
R1="$DNA/01_fastp/SRR19020591.clean.R1.fastq.gz"
R2="$DNA/01_fastp/SRR19020591.clean.R2.fastq.gz"
OUT="$DNA/02b_whole_assembly_DNA_mapping"
THREADS=${THREADS:-16}

mkdir -p "$OUT"
while [[ ! -s "$ASSEMBLY" || ! -s "$R1" || ! -s "$R2" ]]; do
  sleep 60
done

for tool in bowtie2-build bowtie2 samtools; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 2; }
done

if [[ ! -s "$OUT/index.1.bt2" && ! -s "$OUT/index.1.bt2l" ]]; then
  bowtie2-build --threads "$THREADS" "$ASSEMBLY" "$OUT/index" \
    > "$OUT/bowtie2_build.stdout.log" 2> "$OUT/bowtie2_build.stderr.log"
fi

bowtie2 --very-sensitive -x "$OUT/index" -1 "$R1" -2 "$R2" -p "$THREADS" \
  2> "$OUT/SRR19020591.whole_assembly.bowtie2.log" |
  samtools view -@ "$THREADS" -b - |
  samtools sort -@ "$THREADS" -m 2G -o "$OUT/SRR19020591.whole_assembly.sorted.bam" -

samtools index -@ "$THREADS" "$OUT/SRR19020591.whole_assembly.sorted.bam"
samtools flagstat -@ "$THREADS" "$OUT/SRR19020591.whole_assembly.sorted.bam" \
  > "$OUT/SRR19020591.whole_assembly.flagstat.txt"
samtools idxstats "$OUT/SRR19020591.whole_assembly.sorted.bam" \
  > "$OUT/SRR19020591.whole_assembly.idxstats.tsv"
samtools coverage "$OUT/SRR19020591.whole_assembly.sorted.bam" \
  > "$OUT/SRR19020591.whole_assembly.coverage.tsv"

sha256sum "$ASSEMBLY" "$OUT/SRR19020591.whole_assembly.sorted.bam" \
  "$OUT/SRR19020591.whole_assembly.sorted.bam.bai" > "$OUT/SHA256SUMS.txt"
date '+JL DNA whole-assembly mapping complete: %F %T %Z'
