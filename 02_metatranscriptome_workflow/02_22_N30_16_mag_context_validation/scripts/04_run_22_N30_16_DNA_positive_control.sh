#!/usr/bin/env bash
set -euo pipefail

: "${CLEAN_R1:?Set CLEAN_R1 to the audited SRR35213524 clean DNA R1}"
: "${CLEAN_R2:?Set CLEAN_R2 to the audited SRR35213524 clean DNA R2}"
: "${PERSONALIZED_DIR:?Set PERSONALIZED_DIR to finalized reference directory}"
: "${PRODIGAL_GFF:?Set PRODIGAL_GFF to the audited 22_N30_16 Prodigal GFF}"
: "${WHOLE_ASSEMBLY:?Set WHOLE_ASSEMBLY to the audited MEGAHIT final.contigs.fa}"

THREADS=${THREADS:-16}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORK=${WORK:-$(dirname "$PERSONALIZED_DIR")}
MAP="$WORK/05_dna_presence_control_mapping"
DEPTH="$WORK/06_per_base_depth"
FIG="$WORK/07_figures"
IGV="$WORK/08_igv"
TABLE="$WORK/09_final_tables"
MAGGIE=${MAGGIE:-/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra/maggie_expression_validation_20260812}
mkdir -p "$MAP" "$DEPTH" "$FIG" "$IGV" "$TABLE"

for tool in bowtie2 bowtie2-build samtools python3 gzip; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 2; }
done
for file in "$CLEAN_R1" "$CLEAN_R2"; do
  test -s "$file" || { echo "Missing/empty clean DNA read file: $file" >&2; exit 3; }
  gzip -t "$file"
done

OWN="$PERSONALIZED_DIR/22_N30_16_all_K08356_candidates_CDS.fna"
EXTERNAL="$MAGGIE/competitive_DMSOR_plus_treeRefs_CDS.fna"
COMP="$PERSONALIZED_DIR/22_N30_16_personalized_DMSOR_competitive_CDS.fna"
test -s "$OWN" && test -s "$EXTERNAL"
cat "$OWN" "$EXTERNAL" > "$COMP"
WHOLE_REF="$PERSONALIZED_DIR/22_N30_16_whole_assembly_unique_headers.fna"
WHOLE_CROSSWALK="$PERSONALIZED_DIR/22_N30_16_whole_assembly_header_crosswalk.tsv"
python3 "$SCRIPT_DIR/prepare_whole_assembly_reference.py" --assembly "$WHOLE_ASSEMBLY" \
  --out-fasta "$WHOLE_REF" --crosswalk "$WHOLE_CROSSWALK"

declare -A REFERENCES=(
  [whole_assembly_competitive]="$WHOLE_REF"
  [candidate_CDS_competitive]="$COMP"
  [strict_cluster_individual_CDS]="$PERSONALIZED_DIR/22_N30_16_strict_cluster_individual_CDS.fna"
  [strict_cluster_native_regions]="$PERSONALIZED_DIR/22_N30_16_strict_A_B_P-like-1_P-like-2_native_regions.fna"
  [candidate_contigs]="$PERSONALIZED_DIR/22_N30_16_strict_candidate_contigs.fna"
  [all_bins_competitive]="$PERSONALIZED_DIR/22_N30_16_all_bins_competitive.fna"
)

for mode in whole_assembly_competitive all_bins_competitive candidate_CDS_competitive strict_cluster_individual_CDS strict_cluster_native_regions candidate_contigs; do
  ref=${REFERENCES[$mode]}; test -s "$ref" || { echo "Missing reference: $ref" >&2; exit 4; }
  mode_dir="$MAP/$mode"; mkdir -p "$mode_dir"
  index="$mode_dir/reference"; bam="$mode_dir/22_N30_16.${mode}.sorted.bam"
  if [[ ! -e "${index}.1.bt2" && ! -e "${index}.1.bt2l" ]]; then bowtie2-build "$ref" "$index" > "$mode_dir/bowtie2-build.log" 2>&1; fi
  if [[ ! -s "$bam" ]]; then
    bowtie2 -x "$index" -1 "$CLEAN_R1" -2 "$CLEAN_R2" --very-sensitive-local \
      --no-mixed --no-discordant -k 10 -p "$THREADS" 2> "$mode_dir/bowtie2.stderr.log" |
      samtools view -@ "$THREADS" -bS - |
      samtools sort -@ "$THREADS" -o "$bam" -
  fi
  samtools index "$bam"; samtools flagstat "$bam" > "$mode_dir/flagstat.txt"
  samtools idxstats "$bam" > "$mode_dir/idxstats.tsv"
  primary="$mode_dir/22_N30_16.${mode}.primary.sorted.bam"
  mapq20="$mode_dir/22_N30_16.${mode}.primary.MAPQ20.sorted.bam"
  nodup="$mode_dir/22_N30_16.${mode}.primary.no_duplicates.sorted.bam"
  mapq30="$mode_dir/22_N30_16.${mode}.primary.MAPQ30.sorted.bam"
  samtools view -@ "$THREADS" -b -F 2304 "$bam" | samtools sort -@ "$THREADS" -o "$primary" -
  samtools view -@ "$THREADS" -b -F 2304 -q 20 "$bam" | samtools sort -@ "$THREADS" -o "$mapq20" -
  samtools view -@ "$THREADS" -b -F 3328 "$bam" | samtools sort -@ "$THREADS" -o "$nodup" -
  samtools view -@ "$THREADS" -b -F 2304 -q 30 "$bam" | samtools sort -@ "$THREADS" -o "$mapq30" -
  for filtered in "$primary" "$mapq20" "$nodup" "$mapq30"; do samtools index "$filtered"; done
  samtools flagstat "$mapq20" > "$mode_dir/primary.MAPQ20.flagstat.txt"
  samtools idxstats "$mapq20" > "$mode_dir/primary.MAPQ20.idxstats.tsv"
  if [[ "$mode" == "whole_assembly_competitive" ]]; then
    python3 "$SCRIPT_DIR/extract_candidate_contig_depth.py" --bam "$bam" --crosswalk "$WHOLE_CROSSWALK" \
      --strict-genes "$PERSONALIZED_DIR/22_N30_16_strict_cluster_genes.tsv" \
      --out "$DEPTH/22_N30_16.${mode}.candidate_contigs.per_base_depth.tsv.gz"
    python3 "$SCRIPT_DIR/extract_candidate_contig_depth.py" --bam "$mapq20" --crosswalk "$WHOLE_CROSSWALK" \
      --strict-genes "$PERSONALIZED_DIR/22_N30_16_strict_cluster_genes.tsv" \
      --out "$DEPTH/22_N30_16.${mode}.primary.MAPQ20.candidate_contigs.per_base_depth.tsv.gz"
    RAW_DEPTH="$DEPTH/22_N30_16.${mode}.candidate_contigs.per_base_depth.tsv.gz"
    MAIN_DEPTH="$DEPTH/22_N30_16.${mode}.primary.MAPQ20.candidate_contigs.per_base_depth.tsv.gz"
  else
    samtools depth -aa "$bam" | gzip -c > "$DEPTH/22_N30_16.${mode}.per_base_depth.tsv.gz"
    samtools depth -aa "$mapq20" | gzip -c > "$DEPTH/22_N30_16.${mode}.primary.MAPQ20.per_base_depth.tsv.gz"
    RAW_DEPTH="$DEPTH/22_N30_16.${mode}.per_base_depth.tsv.gz"
    MAIN_DEPTH="$DEPTH/22_N30_16.${mode}.primary.MAPQ20.per_base_depth.tsv.gz"
  fi
  python3 "$SCRIPT_DIR/summarize_per_base_depth.py" --depth "$RAW_DEPTH" \
    --mode "$mode" --out "$TABLE/22_N30_16.${mode}.coverage_summary.tsv"
  python3 "$SCRIPT_DIR/summarize_bam_alignment_metrics.py" --bam "$bam" --mode "$mode" \
    --out "$TABLE/22_N30_16.${mode}.alignment_metrics.tsv"
  python3 "$SCRIPT_DIR/summarize_reference_alignment_metrics.py" --bam "$bam" --mode "$mode" \
    --out "$TABLE/22_N30_16.${mode}.per_reference_alignment_metrics.tsv"
  python3 "$SCRIPT_DIR/summarize_per_base_depth.py" --depth "$MAIN_DEPTH" \
    --mode "${mode}_primary_MAPQ20" --out "$TABLE/22_N30_16.${mode}.primary.MAPQ20.coverage_summary.tsv"
  for label_and_bam in "primary:$primary" "primary_no_duplicates:$nodup" "primary_MAPQ30:$mapq30"; do
    label=${label_and_bam%%:*}; filtered=${label_and_bam#*:}
    if [[ "$mode" == "whole_assembly_competitive" ]]; then
      TMP_DEPTH="$DEPTH/22_N30_16.${mode}.${label}.candidate_contigs.per_base_depth.tsv.gz"
      python3 "$SCRIPT_DIR/extract_candidate_contig_depth.py" --bam "$filtered" --crosswalk "$WHOLE_CROSSWALK" \
        --strict-genes "$PERSONALIZED_DIR/22_N30_16_strict_cluster_genes.tsv" --out "$TMP_DEPTH"
      python3 "$SCRIPT_DIR/summarize_per_base_depth.py" --depth "$TMP_DEPTH" \
        --mode "${mode}_${label}" --out "$TABLE/22_N30_16.${mode}.${label}.coverage_summary.tsv"
    else
      samtools depth -aa "$filtered" | python3 "$SCRIPT_DIR/summarize_per_base_depth.py" --depth - \
        --mode "${mode}_${label}" --out "$TABLE/22_N30_16.${mode}.${label}.coverage_summary.tsv"
    fi
    python3 "$SCRIPT_DIR/summarize_reference_alignment_metrics.py" --bam "$filtered" --mode "${mode}_${label}" \
      --out "$TABLE/22_N30_16.${mode}.${label}.per_reference_alignment_metrics.tsv"
  done
  if [[ "$mode" != "all_bins_competitive" ]]; then
    python3 "$SCRIPT_DIR/plot_per_base_depth.py" --depth "$MAIN_DEPTH" \
      --outdir "$FIG/$mode" --title-prefix "DNA presence control (primary MAPQ>=20; not transcriptional evidence)"
  fi
  cp -p "$bam" "$bam.bai" "$primary" "$primary.bai" "$mapq20" "$mapq20.bai" "$ref" "$IGV/"
done

for mode in candidate_CDS_competitive strict_cluster_individual_CDS; do
  python3 "$SCRIPT_DIR/split_per_base_depth.py" \
    --depth "$DEPTH/22_N30_16.${mode}.primary.MAPQ20.per_base_depth.tsv.gz" \
    --outdir "$DEPTH/${mode}_MAPQ20_by_reference"
done

python3 "$SCRIPT_DIR/summarize_MAG_and_CDS_depth.py" \
  --depth "$DEPTH/22_N30_16.all_bins_competitive.primary.MAPQ20.per_base_depth.tsv.gz" \
  --gff "$PRODIGAL_GFF" \
  --bin-manifest "$PERSONALIZED_DIR/22_N30_16_all_bins_reference_manifest.tsv" \
  --strict-genes "$PERSONALIZED_DIR/22_N30_16_strict_cluster_genes.tsv" \
  --out-mag "$TABLE/22_N30_16_host_MAG_DNA_coverage.tsv" \
  --out-cds "$TABLE/22_N30_16_all_binned_CDS_DNA_depth.tsv" \
  --out-strict-cluster "$TABLE/22_N30_16_strict_cluster_genes_all_bins_competitive_DNA_coverage.tsv" \
  --out-candidate-percentile "$TABLE/22_N30_16_IdrA_depth_percentile_within_MAG.tsv"

python3 "$SCRIPT_DIR/count_cluster_boundary_support.py" \
  --bam "$MAP/whole_assembly_competitive/22_N30_16.whole_assembly_competitive.sorted.bam" \
  --boundaries "$PERSONALIZED_DIR/22_N30_16_strict_cluster_boundaries.tsv" \
  --crosswalk "$WHOLE_CROSSWALK" \
  --out "$TABLE/22_N30_16_DNA_cluster_boundary_support.tsv"
python3 "$SCRIPT_DIR/count_cluster_boundary_support.py" \
  --bam "$MAP/whole_assembly_competitive/22_N30_16.whole_assembly_competitive.primary.MAPQ20.sorted.bam" \
  --boundaries "$PERSONALIZED_DIR/22_N30_16_strict_cluster_boundaries.tsv" \
  --crosswalk "$WHOLE_CROSSWALK" \
  --out "$TABLE/22_N30_16_DNA_cluster_boundary_support.MAPQ20.tsv"

python3 "$SCRIPT_DIR/plot_IdrA_HMM_coverage.py" \
  --depth "$DEPTH/22_N30_16.strict_cluster_individual_CDS.primary.MAPQ20.per_base_depth.tsv.gz" \
  --review "$PERSONALIZED_DIR/22_N30_16_final_candidate_classification.tsv" \
  --outdir "$FIG/IdrA_HMM_annotated"

python3 "$SCRIPT_DIR/extract_native_gene_coverage.py" \
  --depth "$DEPTH/22_N30_16.whole_assembly_competitive.primary.MAPQ20.candidate_contigs.per_base_depth.tsv.gz" \
  --crosswalk "$WHOLE_CROSSWALK" \
  --strict-genes "$PERSONALIZED_DIR/22_N30_16_strict_cluster_genes.tsv" \
  --out "$TABLE/22_N30_16_four_gene_native_coverage_primary.tsv"

cat > "$TABLE/22_N30_16_RNA_PENDING.md" <<'EOF'
# 22_N30_16 RNA stage pending

Supplementary Table 1 reports metatranscriptomic sequencing for this physical
sample, but a reliable public RNA run has not been identified. `SRR35213524` is
DNA WGS/METAGENOMIC and was used only for DNA presence-control mapping.

When genuine RNA FASTQ files become available, map them competitively to the
whole-assembly contigs as the primary competitive reference, all bins for host
context, four separate IdrA/IdrB/P-like-1/P-like-2 CDS as an auxiliary view,
and the native cluster region for IGV. Report each gene separately, retain per-base
zero-depth positions, inspect boundary-spanning reads/pairs, and run competitive
DIAMOND blastx. Neighbor-gene support cannot substitute for reliable IdrA CDS
coverage.
EOF

cat > "$TABLE/22_N30_16_DNA_PRESENCE_CONTROL_REPORT.md" <<EOF
# 22_N30_16 DNA presence-control report

This stage uses SRR35213524 metagenomic DNA reads. It is not transcriptional
evidence and cannot establish expression, co-transcription, or operon activity.

## Evidence hierarchy

1. Whole-assembly competitive mapping is the primary candidate-localization
   evidence. Native four-gene coverage is reported in
   \`22_N30_16_four_gene_native_coverage_primary.tsv\`.
2. All-bin competitive mapping describes host-MAG genomic background.
3. Four individual CDS provide auxiliary counts and sensitivity checks only.
4. Native cluster/contig BAMs and junction tables support coordinate and IGV
   inspection. DNA junction evidence is not co-transcription evidence.

## Questions this stage can answer

- Is the candidate present in this DNA sample?
- Does source DNA support the candidate CDS across multiple regions?
- Are the candidate contig and host-MAG assignments supported?
- Is the personalized reference suitable for future RNA analysis?

The corresponding RNA run remains unresolved. No RNA analysis was executed.
EOF

echo "DNA presence-control workflow complete. These outputs are not RNA expression evidence."
