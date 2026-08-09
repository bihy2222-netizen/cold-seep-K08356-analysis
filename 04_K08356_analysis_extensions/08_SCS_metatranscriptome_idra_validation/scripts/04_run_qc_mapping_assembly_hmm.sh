#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:?Usage: bash 03_run_qc_mapping_assembly_hmm.sh MANIFEST.tsv}"
ROOT="${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra}"
RAW_ROOT="${RAW_ROOT:-$ROOT/02_raw_fastq/PRJNA739005_QDN_RNA}"
THREADS="${THREADS:-16}"
ASSEMBLY_THREADS="${ASSEMBLY_THREADS:-24}"
IDRA_HMM="${IDRA_HMM:-}"
STRICT_INDEX="${STRICT_INDEX:-$ROOT/05_idra_reference/idra_strict_20}"
EXTENDED_INDEX="${EXTENDED_INDEX:-$ROOT/05_idra_reference/idra_extended}"
COMPETITIVE_INDEX="${COMPETITIVE_INDEX:-$ROOT/05_idra_reference/DMSOR_competitive}"

mkdir -p "$ROOT"/{03_fastp,04_non_rRNA,06_bowtie2,07_assembly,08_prodigal,09_hmmsearch,10_expression,11_summary,logs,tmp}

FASTP_SUMMARY="$ROOT/11_summary/fastp_summary.tsv"
RRNA_SUMMARY="$ROOT/11_summary/rRNA_removal_summary.tsv"
COMMANDS="$ROOT/11_summary/commands_executed.sh"

[ -s "$FASTP_SUMMARY" ] || printf "sample\trun\tlayout\traw_fastq_R1\traw_fastq_R2\tclean_R1\tclean_R2\tstatus\n" > "$FASTP_SUMMARY"
[ -s "$RRNA_SUMMARY" ] || printf "sample\trun\tclean_R1\tclean_R2\tnon_rRNA_R1\tnon_rRNA_R2\tmethod\tstatus\n" > "$RRNA_SUMMARY"
[ -s "$COMMANDS" ] || printf "# Commands executed for SCS metatranscriptome IdrA validation\n" > "$COMMANDS"

find_fastq_pair() {
  local run="$1"
  local r1 r2
  r1="$(find "$RAW_ROOT" -type f \( -name "${run}_1.fastq.gz" -o -name "${run}_1.fq.gz" -o -name "*${run}*_1.fastq.gz" \) | head -n 1)"
  r2="$(find "$RAW_ROOT" -type f \( -name "${run}_2.fastq.gz" -o -name "${run}_2.fq.gz" -o -name "*${run}*_2.fastq.gz" \) | head -n 1)"
  printf "%s\t%s\n" "$r1" "$r2"
}

run_fastp() {
  local sample="$1" run="$2" r1="$3" r2="$4"
  local outdir="$ROOT/03_fastp/$sample"
  mkdir -p "$outdir"
  local clean1="$outdir/${run}.clean.R1.fastq.gz"
  local clean2="$outdir/${run}.clean.R2.fastq.gz"
  if [ -s "$clean1" ] && [ -s "$clean2" ]; then
    printf "%s\t%s\tPAIRED\t%s\t%s\t%s\t%s\texisting\n" "$sample" "$run" "$r1" "$r2" "$clean1" "$clean2" >> "$FASTP_SUMMARY"
    return
  fi
  echo "fastp $run" >> "$COMMANDS"
  fastp \
    -i "$r1" -I "$r2" \
    -o "$clean1" -O "$clean2" \
    --detect_adapter_for_pe \
    --cut_front --cut_tail \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --n_base_limit 5 \
    --length_required 50 \
    --thread "$THREADS" \
    --html "$outdir/${run}.fastp.html" \
    --json "$outdir/${run}.fastp.json" \
    > "$ROOT/logs/${run}.fastp.log" 2>&1
  printf "%s\t%s\tPAIRED\t%s\t%s\t%s\t%s\tPASS\n" "$sample" "$run" "$r1" "$r2" "$clean1" "$clean2" >> "$FASTP_SUMMARY"
}

run_rrna_removal() {
  local sample="$1" run="$2"
  local clean1="$ROOT/03_fastp/$sample/${run}.clean.R1.fastq.gz"
  local clean2="$ROOT/03_fastp/$sample/${run}.clean.R2.fastq.gz"
  local outdir="$ROOT/04_non_rRNA/$sample"
  mkdir -p "$outdir"
  local nr1="$outdir/${run}.non_rRNA.R1.fastq.gz"
  local nr2="$outdir/${run}.non_rRNA.R2.fastq.gz"
  if [ -s "$nr1" ] && [ -s "$nr2" ]; then
    printf "%s\t%s\t%s\t%s\t%s\t%s\texisting\texisting\n" "$sample" "$run" "$clean1" "$clean2" "$nr1" "$nr2" >> "$RRNA_SUMMARY"
    return
  fi
  if command -v sortmerna >/dev/null 2>&1 && [ -n "${SORTMERNA_REF:-}" ]; then
    echo "sortmerna $run" >> "$COMMANDS"
    sortmerna \
      --ref "$SORTMERNA_REF" \
      --reads "$clean1" --reads "$clean2" \
      --paired_in \
      --fastx \
      --other "$outdir/${run}.non_rRNA" \
      --aligned "$outdir/${run}.rRNA" \
      --workdir "$ROOT/tmp/${run}.sortmerna" \
      --threads "$THREADS" \
      > "$ROOT/logs/${run}.sortmerna.log" 2>&1
    echo "SortMeRNA output naming varies by version; verify and rename non-rRNA paired files if needed." >&2
    printf "%s\t%s\t%s\t%s\t%s\t%s\tsortmerna\tCHECK_OUTPUT_NAMES\n" "$sample" "$run" "$clean1" "$clean2" "$nr1" "$nr2" >> "$RRNA_SUMMARY"
  else
    ln -sf "$clean1" "$nr1"
    ln -sf "$clean2" "$nr2"
    printf "%s\t%s\t%s\t%s\t%s\t%s\tbypass_no_sortmerna_ref\tBYPASS\n" "$sample" "$run" "$clean1" "$clean2" "$nr1" "$nr2" >> "$RRNA_SUMMARY"
  fi
}

run_bowtie_expression() {
  local sample="$1" run="$2" ref_name="$3" index_prefix="$4"
  [ -e "${index_prefix}.1.bt2" ] || [ -e "${index_prefix}.1.bt2l" ] || return 0
  local r1="$ROOT/04_non_rRNA/$sample/${run}.non_rRNA.R1.fastq.gz"
  local r2="$ROOT/04_non_rRNA/$sample/${run}.non_rRNA.R2.fastq.gz"
  local bam="$ROOT/06_bowtie2/${run}.${ref_name}.sorted.bam"
  local mapq="$ROOT/06_bowtie2/${run}.${ref_name}.MAPQ20.bam"
  if [ ! -s "$bam" ]; then
    echo "bowtie2 $run $ref_name" >> "$COMMANDS"
    bowtie2 \
      -x "$index_prefix" \
      -1 "$r1" -2 "$r2" \
      --very-sensitive-local \
      --no-mixed --no-discordant \
      -k 10 \
      -p "$THREADS" \
      2> "$ROOT/logs/${run}.${ref_name}.bowtie2.log" \
      | samtools view -@ 8 -bS - \
      | samtools sort -@ 8 -o "$bam"
    samtools index "$bam"
  fi
  if [ ! -s "$mapq" ]; then
    samtools view -@ 8 -b -q 20 "$bam" > "$mapq"
    samtools index "$mapq"
  fi
  samtools idxstats "$bam" > "$ROOT/10_expression/${run}.${ref_name}.idxstats.tsv"
  samtools depth -aa "$bam" > "$ROOT/10_expression/${run}.${ref_name}.depth.tsv"
  samtools depth -aa "$mapq" > "$ROOT/10_expression/${run}.${ref_name}.MAPQ20.depth.tsv"
}

run_assembly_hmm() {
  local sample="$1" run="$2"
  local r1="$ROOT/04_non_rRNA/$sample/${run}.non_rRNA.R1.fastq.gz"
  local r2="$ROOT/04_non_rRNA/$sample/${run}.non_rRNA.R2.fastq.gz"
  local assembly="$ROOT/07_assembly/$run"
  if [ ! -s "$assembly/final.contigs.fa" ]; then
    if [ -d "$assembly" ]; then
      echo "Assembly directory exists without final.contigs.fa: $assembly" >&2
      return 1
    fi
    echo "megahit $run" >> "$COMMANDS"
    megahit \
      -1 "$r1" -2 "$r2" \
      --presets meta-sensitive \
      --min-contig-len 300 \
      -t "$ASSEMBLY_THREADS" \
      -m 0.7 \
      -o "$assembly" \
      > "$ROOT/logs/${run}.megahit.log" 2>&1
  fi
  if [ ! -s "$ROOT/08_prodigal/${run}.proteins.faa" ]; then
    echo "prodigal $run" >> "$COMMANDS"
    prodigal \
      -i "$assembly/final.contigs.fa" \
      -a "$ROOT/08_prodigal/${run}.proteins.faa" \
      -d "$ROOT/08_prodigal/${run}.CDS.fna" \
      -f gff \
      -o "$ROOT/08_prodigal/${run}.genes.gff" \
      -p meta \
      > "$ROOT/logs/${run}.prodigal.log" 2>&1
  fi
  if [ -n "$IDRA_HMM" ] && [ -s "$IDRA_HMM" ]; then
    echo "hmmsearch $run" >> "$COMMANDS"
    hmmsearch \
      --cpu "$THREADS" \
      --tblout "$ROOT/09_hmmsearch/${run}.idra.tbl" \
      --domtblout "$ROOT/09_hmmsearch/${run}.idra.domtbl" \
      "$IDRA_HMM" \
      "$ROOT/08_prodigal/${run}.proteins.faa" \
      > "$ROOT/09_hmmsearch/${run}.idra.hmmsearch.txt"
  else
    echo "IDRA_HMM not set or missing; skipped HMM search for $run" | tee -a "$ROOT/logs/${run}.hmmsearch.SKIPPED.log"
  fi
}

python3 - "$MANIFEST" <<'PY' > "$ROOT/tmp/manifest_samples.tsv"
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
for r in rows:
    print("\t".join([r.get("Run",""), r.get("SampleName","sample"), r.get("LibraryLayout","")]))
PY

while IFS=$'\t' read -r run sample layout; do
  [ -n "$run" ] || continue
  pair="$(find_fastq_pair "$run")"
  r1="$(printf "%s" "$pair" | cut -f1)"
  r2="$(printf "%s" "$pair" | cut -f2)"
  if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
    echo "Missing paired FASTQ for $run ($sample); skipping for now." | tee -a "$ROOT/logs/missing_fastq.log"
    continue
  fi
  gzip -t "$r1"
  gzip -t "$r2"
  run_fastp "$sample" "$run" "$r1" "$r2"
  run_rrna_removal "$sample" "$run"
  run_bowtie_expression "$sample" "$run" "strict" "$STRICT_INDEX"
  run_bowtie_expression "$sample" "$run" "extended" "$EXTENDED_INDEX"
  run_bowtie_expression "$sample" "$run" "competitive" "$COMPETITIVE_INDEX"
  run_assembly_hmm "$sample" "$run"
done < "$ROOT/tmp/manifest_samples.tsv"

echo "QC/mapping/assembly/HMM workflow finished: $(date)"
