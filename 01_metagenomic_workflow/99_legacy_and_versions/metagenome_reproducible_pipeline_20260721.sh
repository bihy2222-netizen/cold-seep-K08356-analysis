#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Reproducible metagenome pipeline for cold-seep habitat comparison
#
# This script is a corrected and standardized replacement for the personal
# command notes in:
#   宏基因组分析流程全部脚本整理.R
#
# The original file is intentionally left unchanged.
#
# Default behavior is DRY_RUN=1: commands are printed, not executed.
# After checking paths and samples.tsv, run selected steps with DRY_RUN=0.
#
# Example:
#   bash metagenome_reproducible_pipeline_20260721.sh init
#   bash metagenome_reproducible_pipeline_20260721.sh check
#   DRY_RUN=0 bash metagenome_reproducible_pipeline_20260721.sh qc
#
# Required samples.tsv columns:
#   sample  habitat  site  depth  r1  r2  srr
#
# Notes:
#   - r1/r2 are cleanable raw FASTQ paths for local data.
#   - srr can be NA if data are already downloaded.
#   - Use tabs, not commas.
###############################################################################

#######################################
# 0. Configuration
#######################################

PROJECT_DIR="${PROJECT_DIR:-$(pwd)/metagenome_project}"
SAMPLES_TSV="${SAMPLES_TSV:-${PROJECT_DIR}/metadata/samples.tsv}"
THREADS="${THREADS:-48}"
MAX_JOBS="${MAX_JOBS:-4}"
DRY_RUN="${DRY_RUN:-1}"
SEED="${SEED:-12345}"

# Conda activation: set CONDA_SH if conda is not available in non-interactive bash.
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"

# Databases and references. Edit these before running real analyses.
KOFAM_EXEC="${KOFAM_EXEC:-/home/ps/Research/databases/other_database/kofamscan/kofam_scan-1.3.0/exec_annotation}"
KOFAM_PROFILES="${KOFAM_PROFILES:-/home/ps/Research/databases/other_database/kofamscan/profiles}"
KOFAM_KO_LIST="${KOFAM_KO_LIST:-/home/ps/Research/databases/other_database/kofamscan/ko_list}"

PFAM00384_HMM="${PFAM00384_HMM:-${PROJECT_DIR}/references/PF00384.hmm}"
AIOA_SEED_FAA="${AIOA_SEED_FAA:-${PROJECT_DIR}/references/aioA_seed.faa}"
IDRA_SEED_FAA="${IDRA_SEED_FAA:-${PROJECT_DIR}/references/idrA_seed.faa}"
DMSOR_REFERENCE_FAA="${DMSOR_REFERENCE_FAA:-${PROJECT_DIR}/references/dmsor_family_references.faa}"

LOG_DIR="${PROJECT_DIR}/logs"

ensure_project_dirs() {
  mkdir -p "${PROJECT_DIR}"/{metadata,logs,00_raw/sra,00_raw/fastq,01_qc/fastqc_raw,01_qc/fastqc_clean,01_qc/fastp,02_clean,03_assembly,03_assembly_qc,04_metaphlan,05_binning,06_genes/{CDS,AA,catalog},07_abundance/{gene,MAG},08_annotation/{kofam,k08356,hmm,neighborhood},09_MAGs/{all_refined_bins,filtered_MAGs,dRep_out,checkm2,gunc,gtdbtk},10_statistics,reports,references}
}

#######################################
# 1. Helpers
#######################################

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

run() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    printf '[DRY_RUN] %s\n' "$*"
  else
    log "RUN: $*"
    eval "$@"
  fi
}

activate_env() {
  local env_name="$1"
  if [[ -f "${CONDA_SH}" ]]; then
    # shellcheck source=/dev/null
    source "${CONDA_SH}"
  fi
  run "conda activate ${env_name}"
}

require_samples() {
  if [[ ! -f "${SAMPLES_TSV}" ]]; then
    cat >&2 <<EOF
Missing samples.tsv:
  ${SAMPLES_TSV}

Create it with columns:
  sample<TAB>habitat<TAB>site<TAB>depth<TAB>r1<TAB>r2<TAB>srr

Use:
  bash $0 init
EOF
    exit 1
  fi
}

sample_rows() {
  require_samples
  awk 'BEGIN{FS=OFS="\t"} NR > 1 && $1 !~ /^#/ {print}' "${SAMPLES_TSV}"
}

sample_col() {
  local col="$1"
  awk -v target="${col}" '
    BEGIN{FS=OFS="\t"}
    NR == 1 {
      for (i = 1; i <= NF; i++) if ($i == target) idx = i
      if (!idx) exit 2
      next
    }
    NR > 1 && $1 !~ /^#/ {print $idx}
  ' "${SAMPLES_TSV}"
}

#######################################
# 2. Initialization and checks
#######################################

init_project() {
  ensure_project_dirs
  mkdir -p "${PROJECT_DIR}/metadata"
  if [[ ! -f "${SAMPLES_TSV}" ]]; then
    cat > "${SAMPLES_TSV}" <<'EOF'
sample	habitat	site	depth	r1	r2	srr
NS_0-6	NS	NS	0-6	/path/to/NS_0-6_R1.fastq.gz	/path/to/NS_0-6_R2.fastq.gz	NA
ES_2_0-6	ES	ES_2	0-6	/path/to/ES_2_0-6_R1.fastq.gz	/path/to/ES_2_0-6_R2.fastq.gz	NA
EOF
  fi

  cat > "${PROJECT_DIR}/metadata/reference_protein_manifest.tsv" <<'EOF'
sequence_id	family	accession	evidence_level	source	note
example_AioA	AioA	NA	experimentally_validated	literature	Replace with curated AioA reference
example_IdrA	IdrA	NA	experimentally_validated	literature	Include Denitromonas IR-12 and Pseudomonas SCT IdrA
example_ArxA	ArxA	NA	curated	literature	Negative/near-neighbor DMSOR family
example_ArrA	ArrA	NA	curated	literature	Negative/near-neighbor DMSOR family
example_DmsA	DmsA	NA	curated	literature	Outgroup DMSOR family
example_NarG	NarG	NA	curated	literature	Outgroup DMSOR family
example_NapA	NapA	NA	curated	literature	Outgroup DMSOR family
example_PsrA	PsrA	NA	curated	literature	Outgroup DMSOR family
EOF

  cat > "${PROJECT_DIR}/metadata/clade_assignment_evidence.tsv" <<'EOF'
sequence_id	habitat	host_MAG	AioA_score	IdrA_score	AioA_cov	IdrA_cov	score_delta	tree_clade	UFboot	neighborhood	contig_edge	final_class	confidence
EOF

  log "Initialized project folders and template metadata under ${PROJECT_DIR}"
}

check_inputs() {
  require_samples
  log "Checking sample table columns"
  awk 'BEGIN{FS="\t"; ok=1}
    NR == 1 {
      required = "sample habitat site depth r1 r2 srr"
      n = split(required, req, " ")
      for (i = 1; i <= NF; i++) seen[$i] = 1
      for (i = 1; i <= n; i++) {
        if (!(req[i] in seen)) {
          print "Missing column: " req[i] > "/dev/stderr"
          ok = 0
        }
      }
      exit ok ? 0 : 1
    }' "${SAMPLES_TSV}"

  log "Checking FASTQ paths listed in samples.tsv"
  awk 'BEGIN{FS="\t"; missing=0}
    NR == 1 {
      for (i = 1; i <= NF; i++) idx[$i] = i
      next
    }
    NR > 1 && $1 !~ /^#/ {
      r1 = $idx["r1"]; r2 = $idx["r2"]
      if (r1 != "NA" && system("[ -f \"" r1 "\" ]") != 0) {
        print "Missing R1 for " $idx["sample"] ": " r1 > "/dev/stderr"; missing=1
      }
      if (r2 != "NA" && system("[ -f \"" r2 "\" ]") != 0) {
        print "Missing R2 for " $idx["sample"] ": " r2 > "/dev/stderr"; missing=1
      }
    }
    END{exit missing ? 1 : 0}' "${SAMPLES_TSV}"

  log "Input check completed"
}

record_versions() {
  local out="${PROJECT_DIR}/metadata/software_versions.tsv"
  run "printf 'software\tversion\n' > '${out}'"
  for cmd in fastqc multiqc fastp fasterq-dump vdb-validate megahit quast.py metaphlan metawrap prodigal cd-hit-est cd-hit salmon checkm2 gunc dRep gtdb-tk trimal iqtree2 iqtree hmmbuild hmmpress hmmsearch mafft seqkit coverm; do
    run "printf '${cmd}\t' >> '${out}'; (${cmd} --version 2>&1 | head -n 1 || true) >> '${out}'"
  done
}

#######################################
# 3. Download and read QC
#######################################

download_sra() {
  activate_env sratools
  require_samples
  sample_rows | while IFS=$'\t' read -r sample habitat site depth r1 r2 srr; do
    [[ "${srr}" == "NA" || -z "${srr}" ]] && continue
    run "prefetch '${srr}' -O '${PROJECT_DIR}/00_raw/sra' > '${LOG_DIR}/${sample}.prefetch.log' 2>&1"
    run "vdb-validate '${PROJECT_DIR}/00_raw/sra/${srr}' > '${LOG_DIR}/${sample}.vdb_validate.log' 2>&1"
    run "fasterq-dump '${PROJECT_DIR}/00_raw/sra/${srr}' --split-files --threads 16 --outdir '${PROJECT_DIR}/00_raw/fastq' > '${LOG_DIR}/${sample}.fasterq_dump.log' 2>&1"
  done
  run "md5sum '${PROJECT_DIR}'/00_raw/fastq/*.fastq* > '${PROJECT_DIR}/metadata/raw_fastq.md5' || true"
}

qc_raw() {
  activate_env fastqc
  run "fastqc '${PROJECT_DIR}'/00_raw/fastq/*.fastq* -o '${PROJECT_DIR}/01_qc/fastqc_raw' -t '${THREADS}' > '${LOG_DIR}/fastqc_raw.log' 2>&1"
  activate_env multiqc
  run "multiqc '${PROJECT_DIR}/01_qc/fastqc_raw' -o '${PROJECT_DIR}/01_qc' > '${LOG_DIR}/multiqc_raw.log' 2>&1"
}

fastp_clean() {
  activate_env fastp
  require_samples
  sample_rows | while IFS=$'\t' read -r sample habitat site depth r1 r2 srr; do
    [[ "${r1}" == "NA" || "${r2}" == "NA" ]] && continue
    run "fastp -i '${r1}' -I '${r2}' -o '${PROJECT_DIR}/02_clean/${sample}_R1.clean.fastq.gz' -O '${PROJECT_DIR}/02_clean/${sample}_R2.clean.fastq.gz' --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50 --thread 16 --html '${PROJECT_DIR}/01_qc/fastp/${sample}.fastp.html' --json '${PROJECT_DIR}/01_qc/fastp/${sample}.fastp.json' > '${LOG_DIR}/${sample}.fastp.log' 2>&1"
  done
  run "md5sum '${PROJECT_DIR}'/02_clean/*.clean.fastq.gz > '${PROJECT_DIR}/metadata/clean_fastq.md5' || true"
}

qc_clean() {
  activate_env fastqc
  run "fastqc '${PROJECT_DIR}'/02_clean/*.clean.fastq.gz -o '${PROJECT_DIR}/01_qc/fastqc_clean' -t '${THREADS}' > '${LOG_DIR}/fastqc_clean.log' 2>&1"
  activate_env multiqc
  run "multiqc '${PROJECT_DIR}/01_qc/fastqc_clean' '${PROJECT_DIR}/01_qc/fastp' -o '${PROJECT_DIR}/01_qc' > '${LOG_DIR}/multiqc_clean.log' 2>&1"
}

#######################################
# 4. Assembly and community composition
#######################################

assembly_megahit() {
  activate_env megahit
  sample_col sample | while read -r sample; do
    [[ -z "${sample}" ]] && continue
    run "megahit -1 '${PROJECT_DIR}/02_clean/${sample}_R1.clean.fastq.gz' -2 '${PROJECT_DIR}/02_clean/${sample}_R2.clean.fastq.gz' --k-min 21 --k-max 141 --k-step 10 --min-contig-len 1000 -t '${THREADS}' -o '${PROJECT_DIR}/03_assembly/${sample}' > '${LOG_DIR}/${sample}.megahit.log' 2>&1"
  done
}

assembly_quast() {
  activate_env quast
  sample_col sample | while read -r sample; do
    run "quast.py '${PROJECT_DIR}/03_assembly/${sample}/final.contigs.fa' -o '${PROJECT_DIR}/03_assembly_qc/${sample}' -t 16 > '${LOG_DIR}/${sample}.quast.log' 2>&1"
  done
}

metaphlan_profile() {
  activate_env metaphlan
  sample_col sample | while read -r sample; do
    run "metaphlan '${PROJECT_DIR}/02_clean/${sample}_R1.clean.fastq.gz','${PROJECT_DIR}/02_clean/${sample}_R2.clean.fastq.gz' --input_type fastq --nproc 32 --bowtie2out '${PROJECT_DIR}/04_metaphlan/${sample}.bowtie2.bz2' -o '${PROJECT_DIR}/04_metaphlan/${sample}.metaphlan.tsv' > '${LOG_DIR}/${sample}.metaphlan.log' 2>&1"
  done
  run "cd '${PROJECT_DIR}/04_metaphlan' && merge_metaphlan_tables.py *.metaphlan.tsv > merged_abundance_table.txt"
}

#######################################
# 5. Binning, MAG quality, dRep, GTDB-Tk
#######################################

metawrap_binning() {
  activate_env metawrap
  sample_col sample | while read -r sample; do
    run "metawrap binning -o '${PROJECT_DIR}/05_binning/${sample}_bins' -t '${THREADS}' -a '${PROJECT_DIR}/03_assembly/${sample}/final.contigs.fa' --metabat2 --maxbin2 --concoct '${PROJECT_DIR}/02_clean/${sample}_R1.clean.fastq.gz' '${PROJECT_DIR}/02_clean/${sample}_R2.clean.fastq.gz' > '${LOG_DIR}/${sample}.metawrap_binning.log' 2>&1"
  done
}

metawrap_refinement() {
  activate_env metawrap
  sample_col sample | while read -r sample; do
    run "metawrap bin_refinement -o '${PROJECT_DIR}/05_binning/${sample}_bins/BIN_REFINEMENT' -t '${THREADS}' -A '${PROJECT_DIR}/05_binning/${sample}_bins/metabat2_bins' -B '${PROJECT_DIR}/05_binning/${sample}_bins/maxbin2_bins' -C '${PROJECT_DIR}/05_binning/${sample}_bins/concoct_bins' -c 50 -x 10 > '${LOG_DIR}/${sample}.metawrap_refinement.log' 2>&1"
  done
}

combine_bin_stats() {
  local out="${PROJECT_DIR}/09_MAGs/all_combined_bins.stats"
  run "printf 'sample\tbin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner\n' > '${out}'"
  run "find '${PROJECT_DIR}/05_binning' -path '*/BIN_REFINEMENT/metawrap_50_10_bins.stats' -print | while read -r file; do sample_dir=\$(basename \"\$(dirname \"\$(dirname \"\$file\")\")\"); sample_name=\${sample_dir%_bins}; awk -v sample=\"\$sample_name\" 'NR > 1 {gsub(/^bin\\./, \"\", \$1); printf \"%s\\t%s\", sample, \$1; for (i = 2; i <= NF; i++) printf \"\\t%s\", \$i; printf \"\\n\"}' OFS='\\t' \"\$file\" >> '${out}'; done"
}

collect_refined_bins() {
  # Dry-run this first. It copies bins and prefixes sample names to avoid ID collisions.
  sample_col sample | while read -r sample; do
    run "find '${PROJECT_DIR}/05_binning/${sample}_bins/BIN_REFINEMENT' -path '*metawrap_50_10_bins/*.fa' -print | while read -r bin; do base=\$(basename \"\$bin\" .fa); cp \"\$bin\" '${PROJECT_DIR}/09_MAGs/all_refined_bins/${sample}_'\"\$base\"'.fasta'; done"
  done
}

mag_quality_filter() {
  activate_env checkm2
  run "checkm2 predict --input '${PROJECT_DIR}/09_MAGs/all_refined_bins' --output-directory '${PROJECT_DIR}/09_MAGs/checkm2' --threads '${THREADS}' --extension fasta > '${LOG_DIR}/checkm2.log' 2>&1"

  activate_env gunc
  run "gunc run --input_dir '${PROJECT_DIR}/09_MAGs/all_refined_bins' --file_suffix .fasta --out_dir '${PROJECT_DIR}/09_MAGs/gunc' --threads 32 > '${LOG_DIR}/gunc.log' 2>&1"

  # Filter with CheckM2 table. Adjust column names if CheckM2 output format differs.
  run "awk 'BEGIN{FS=OFS=\"\\t\"} NR==1 {for(i=1;i<=NF;i++) h[\$i]=i; print; next} NR>1 && \$(h[\"Completeness\"]) >= 50 && \$(h[\"Contamination\"]) <= 10 {print}' '${PROJECT_DIR}/09_MAGs/checkm2/quality_report.tsv' > '${PROJECT_DIR}/09_MAGs/checkm2/quality_report.filtered.tsv'"
  run "awk 'BEGIN{FS=\"\\t\"} NR>1 {print \$1}' '${PROJECT_DIR}/09_MAGs/checkm2/quality_report.filtered.tsv' | while read -r mag; do cp '${PROJECT_DIR}/09_MAGs/all_refined_bins/'\"\$mag\"'.fasta' '${PROJECT_DIR}/09_MAGs/filtered_MAGs/' || true; done"
}

drep_all_samples() {
  # Do not use --ignoreGenomeQuality. Use explicit ANI and quality thresholds.
  activate_env drep
  run "dRep dereplicate '${PROJECT_DIR}/09_MAGs/dRep_out' -g '${PROJECT_DIR}/09_MAGs/filtered_MAGs/'*.fasta -p '${THREADS}' -comp 50 -con 10 -sa 0.95 -nc 0.30 > '${LOG_DIR}/drep.log' 2>&1"
}

gtdbtk_classify_and_tree() {
  activate_env gtdbtk
  run "gtdb-tk classify_wf --genome_dir '${PROJECT_DIR}/09_MAGs/dRep_out/dereplicated_genomes' --out_dir '${PROJECT_DIR}/09_MAGs/gtdbtk/classify_wf_out' --cpus '${THREADS}' --extension fasta > '${LOG_DIR}/gtdbtk_classify.log' 2>&1"

  # Important correction: GTDB marker MSA is already aligned. Do not run MAFFT again.
  activate_env trimal
  run "cd '${PROJECT_DIR}/09_MAGs/gtdbtk/classify_wf_out/align' && gunzip -f -k gtdbtk.*.user_msa.fasta.gz || true"
  run "cd '${PROJECT_DIR}/09_MAGs/gtdbtk/classify_wf_out/align' && if [[ -f gtdbtk.bac120.user_msa.fasta ]]; then trimal -in gtdbtk.bac120.user_msa.fasta -out gtdbtk.bac120.user_msa.trimmed.fasta -gappyout -fasta; fi"
  run "cd '${PROJECT_DIR}/09_MAGs/gtdbtk/classify_wf_out/align' && if [[ -f gtdbtk.ar53.user_msa.fasta ]]; then trimal -in gtdbtk.ar53.user_msa.fasta -out gtdbtk.ar53.user_msa.trimmed.fasta -gappyout -fasta; fi"

  activate_env iqtree
  run "cd '${PROJECT_DIR}/09_MAGs/gtdbtk/classify_wf_out/align' && if command -v iqtree2 >/dev/null 2>&1; then IQTREE=iqtree2; else IQTREE=iqtree; fi; if [[ -f gtdbtk.bac120.user_msa.trimmed.fasta ]]; then \${IQTREE} -s gtdbtk.bac120.user_msa.trimmed.fasta -m MFP -B 1000 --alrt 1000 -T AUTO --seed '${SEED}' --prefix gtdbtk.bac120; fi"
  run "cd '${PROJECT_DIR}/09_MAGs/gtdbtk/classify_wf_out/align' && if command -v iqtree2 >/dev/null 2>&1; then IQTREE=iqtree2; else IQTREE=iqtree; fi; if [[ -f gtdbtk.ar53.user_msa.trimmed.fasta ]]; then \${IQTREE} -s gtdbtk.ar53.user_msa.trimmed.fasta -m MFP -B 1000 --alrt 1000 -T AUTO --seed '${SEED}' --prefix gtdbtk.ar53; fi"
}

#######################################
# 6. Gene prediction, unified catalog, abundance
#######################################

predict_genes() {
  activate_env prodigal
  sample_col sample | while read -r sample; do
    run "prodigal -i '${PROJECT_DIR}/03_assembly/${sample}/final.contigs.fa' -d '${PROJECT_DIR}/06_genes/CDS/${sample}_cds.fa' -a '${PROJECT_DIR}/06_genes/AA/${sample}_AA.faa' -p meta -m > '${LOG_DIR}/${sample}.prodigal.log' 2>&1"
  done
}

build_gene_catalog() {
  # Unified non-redundant catalog for cross-habitat comparison.
  activate_env cd-hit
  run "cat '${PROJECT_DIR}'/06_genes/CDS/*_cds.fa > '${PROJECT_DIR}/06_genes/catalog/all_samples_cds.fa'"
  run "cd-hit-est -i '${PROJECT_DIR}/06_genes/catalog/all_samples_cds.fa' -o '${PROJECT_DIR}/06_genes/catalog/all_samples_cds.nr95.fa' -c 0.95 -aS 0.9 -G 0 -M 0 -g 1 -T '${THREADS}' > '${LOG_DIR}/cdhit_catalog_cds.log' 2>&1"
  run "cat '${PROJECT_DIR}'/06_genes/AA/*_AA.faa > '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.faa'"
  run "cd-hit -i '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.faa' -o '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.nr95.faa' -c 0.95 -aS 0.9 -M 0 -T '${THREADS}' > '${LOG_DIR}/cdhit_catalog_AA.log' 2>&1"
}

salmon_gene_abundance() {
  activate_env salmon
  run "salmon index -t '${PROJECT_DIR}/06_genes/catalog/all_samples_cds.nr95.fa' -i '${PROJECT_DIR}/07_abundance/gene/salmon_catalog_index' -p 24 > '${LOG_DIR}/salmon_catalog_index.log' 2>&1"
  sample_col sample | while read -r sample; do
    run "salmon quant -i '${PROJECT_DIR}/07_abundance/gene/salmon_catalog_index' -l A -1 '${PROJECT_DIR}/02_clean/${sample}_R1.clean.fastq.gz' -2 '${PROJECT_DIR}/02_clean/${sample}_R2.clean.fastq.gz' -o '${PROJECT_DIR}/07_abundance/gene/${sample}' --validateMappings --meta -p '${THREADS}' > '${LOG_DIR}/${sample}.salmon_catalog_quant.log' 2>&1"
  done
}

mag_abundance_coverm() {
  activate_env coverm
  sample_col sample | while read -r sample; do
    run "coverm genome --coupled '${PROJECT_DIR}/02_clean/${sample}_R1.clean.fastq.gz' '${PROJECT_DIR}/02_clean/${sample}_R2.clean.fastq.gz' --genome-fasta-directory '${PROJECT_DIR}/09_MAGs/dRep_out/dereplicated_genomes' --methods relative_abundance trimmed_mean covered_fraction --min-read-percent-identity 95 --min-read-aligned-percent 75 --min-covered-fraction 0.5 --threads '${THREADS}' --output-file '${PROJECT_DIR}/07_abundance/MAG/${sample}.mag_abundance.tsv' > '${LOG_DIR}/${sample}.coverm_MAG.log' 2>&1"
  done
}

#######################################
# 7. Functional annotation and K08356 closed loop
#######################################

kofamscan_annotation() {
  activate_env kofamscan
  run "'${KOFAM_EXEC}' --cpu '${THREADS}' -f detail-tsv --profile '${KOFAM_PROFILES}' --ko-list '${KOFAM_KO_LIST}' --tmp-dir '${PROJECT_DIR}/08_annotation/kofam/tmp' -o '${PROJECT_DIR}/08_annotation/kofam/all_samples_AA.nr95.kofamscan.tsv' '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.nr95.faa' > '${LOG_DIR}/kofamscan.log' 2>&1"
}

k08356_candidate_search() {
  activate_env hmmer
  run "hmmpress '${PFAM00384_HMM}' || true"
  run "hmmsearch --cpu '${THREADS}' --domtblout '${PROJECT_DIR}/08_annotation/k08356/PF00384.domtblout' --tblout '${PROJECT_DIR}/08_annotation/k08356/PF00384.tblout' '${PFAM00384_HMM}' '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.nr95.faa' > '${LOG_DIR}/PF00384.hmmsearch.log' 2>&1"

  # KOfam K08356 candidates, if KOfam output includes KO identifiers.
  run "awk 'BEGIN{FS=OFS=\"\\t\"} /K08356/ {print}' '${PROJECT_DIR}/08_annotation/kofam/all_samples_AA.nr95.kofamscan.tsv' > '${PROJECT_DIR}/08_annotation/k08356/K08356.kofam_hits.tsv' || true"

  # TODO: parse domtblout to extract candidate FASTA IDs. Keep this as an explicit step
  # because HMMER domtblout column layout and ID conventions must be checked.
  run "printf 'candidate_id\\tsource\\tnote\\n' > '${PROJECT_DIR}/08_annotation/k08356/k08356_candidates.tsv'"
}

build_aioa_idra_hmms() {
  activate_env mafft
  run "mafft --localpair --maxiterate 1000 '${IDRA_SEED_FAA}' > '${PROJECT_DIR}/08_annotation/hmm/idrA_seed.aln.faa'"
  run "mafft --localpair --maxiterate 1000 '${AIOA_SEED_FAA}' > '${PROJECT_DIR}/08_annotation/hmm/aioA_seed.aln.faa'"

  activate_env trimal
  run "trimal -in '${PROJECT_DIR}/08_annotation/hmm/idrA_seed.aln.faa' -out '${PROJECT_DIR}/08_annotation/hmm/idrA_seed.trim.faa' -gappyout"
  run "trimal -in '${PROJECT_DIR}/08_annotation/hmm/aioA_seed.aln.faa' -out '${PROJECT_DIR}/08_annotation/hmm/aioA_seed.trim.faa' -gappyout"

  activate_env hmmer
  run "hmmbuild '${PROJECT_DIR}/08_annotation/hmm/idrA.hmm' '${PROJECT_DIR}/08_annotation/hmm/idrA_seed.trim.faa' > '${LOG_DIR}/hmmbuild_idrA.log' 2>&1"
  run "hmmbuild '${PROJECT_DIR}/08_annotation/hmm/aioA.hmm' '${PROJECT_DIR}/08_annotation/hmm/aioA_seed.trim.faa' > '${LOG_DIR}/hmmbuild_aioA.log' 2>&1"
  run "hmmpress '${PROJECT_DIR}/08_annotation/hmm/idrA.hmm' || true"
  run "hmmpress '${PROJECT_DIR}/08_annotation/hmm/aioA.hmm' || true"
}

aioa_idra_dual_hmmsearch() {
  activate_env hmmer
  local candidates="${PROJECT_DIR}/08_annotation/k08356/k08356_candidates.faa"
  run "hmmsearch --cpu '${THREADS}' --domtblout '${PROJECT_DIR}/08_annotation/k08356/idrA.domtblout' --tblout '${PROJECT_DIR}/08_annotation/k08356/idrA.tblout' '${PROJECT_DIR}/08_annotation/hmm/idrA.hmm' '${candidates}' > '${LOG_DIR}/idrA.hmmsearch.log' 2>&1"
  run "hmmsearch --cpu '${THREADS}' --domtblout '${PROJECT_DIR}/08_annotation/k08356/aioA.domtblout' --tblout '${PROJECT_DIR}/08_annotation/k08356/aioA.tblout' '${PROJECT_DIR}/08_annotation/hmm/aioA.hmm' '${candidates}' > '${LOG_DIR}/aioA.hmmsearch.log' 2>&1"
}

k08356_reference_tree() {
  # Include AioA, IdrA, ArxA, ArrA, DmsA/DorA, NarG/NapA, PsrA/PhsA, TtrA, SerA,
  # and other reliable MopB/DMSOR outgroups. Do not classify with only AioA/IdrA.
  activate_env seqkit
  run "seqkit stats '${PROJECT_DIR}/08_annotation/k08356/k08356_candidates.faa' '${DMSOR_REFERENCE_FAA}' > '${PROJECT_DIR}/08_annotation/k08356/tree_input.seqkit_stats.tsv'"
  run "seqkit fx2tab -n -l '${PROJECT_DIR}/08_annotation/k08356/k08356_candidates.faa' > '${PROJECT_DIR}/08_annotation/k08356/candidate_lengths.tsv'"

  activate_env mafft
  run "cat '${PROJECT_DIR}/08_annotation/k08356/k08356_candidates.faa' '${DMSOR_REFERENCE_FAA}' > '${PROJECT_DIR}/08_annotation/k08356/k08356_plus_DMSOR_refs.faa'"
  run "mafft --localpair --maxiterate 1000 '${PROJECT_DIR}/08_annotation/k08356/k08356_plus_DMSOR_refs.faa' > '${PROJECT_DIR}/08_annotation/k08356/k08356_plus_DMSOR_refs.aln.faa' 2> '${LOG_DIR}/k08356_mafft.log'"

  activate_env trimal
  run "trimal -in '${PROJECT_DIR}/08_annotation/k08356/k08356_plus_DMSOR_refs.aln.faa' -out '${PROJECT_DIR}/08_annotation/k08356/k08356_plus_DMSOR_refs.trimmed.faa' -gappyout"

  activate_env iqtree
  run "cd '${PROJECT_DIR}/08_annotation/k08356' && if command -v iqtree2 >/dev/null 2>&1; then IQTREE=iqtree2; else IQTREE=iqtree; fi; \${IQTREE} -s k08356_plus_DMSOR_refs.trimmed.faa -m MFP -B 1000 --alrt 1000 -T AUTO --seed '${SEED}' --prefix k08356_DMSOR_reference_tree"
}

gene_neighborhood_template() {
  # Placeholder extraction table. Use Prodigal GFF/Prokka/DRAM output for real extraction.
  # Record contig edge effects so missing operon genes are not overinterpreted.
  local out="${PROJECT_DIR}/08_annotation/neighborhood/k08356_gene_neighborhood.tsv"
  run "printf 'target_id\\tcontig\\tcontig_length\\ttarget_start\\ttarget_end\\tdistance_to_left_edge\\tdistance_to_right_edge\\tupstream_10_orfs\\tdownstream_10_orfs\\taioB_or_idrB\\tP1_P2\\telectron_transfer_gene\\tarsenic_resistance_gene\\ttransposase_integrase\\tneighborhood_class\\n' > '${out}'"
}

clade_assignment_evidence_template() {
  # Final class must combine HMM, tree and neighborhood evidence.
  # Classes: canonical_AioA_associated, IdrA_associated, AioA_like_associated,
  # unknown_DMSOR. Ambiguous/bootstrap-unstable cases should be unknown.
  local out="${PROJECT_DIR}/metadata/clade_assignment_evidence.tsv"
  run "printf 'sequence_id\\thabitat\\thost_MAG\\tAioA_score\\tIdrA_score\\tAioA_cov\\tIdrA_cov\\tscore_delta\\ttree_clade\\tUFboot\\tneighborhood\\tcontig_edge\\tfinal_class\\tconfidence\\n' > '${out}'"
}

#######################################
# 8. Statistics skeleton
#######################################

write_statistics_template() {
  local rfile="${PROJECT_DIR}/10_statistics/08_statistics_template.R"
  cat > "${rfile}" <<'EOF'
# Statistics template for four-habitat metagenome comparison
# Fill input tables before running.

packages <- c("tidyverse", "lme4", "lmerTest", "emmeans", "effectsize")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)

abund_file <- "../07_abundance/gene/k08356_clade_abundance.tsv"
meta_file <- "../metadata/samples.tsv"

abund <- readr::read_tsv(abund_file, show_col_types = FALSE)
meta <- readr::read_tsv(meta_file, show_col_types = FALSE)

dat <- abund %>%
  left_join(meta, by = "sample") %>%
  mutate(
    habitat = factor(habitat),
    site = factor(site),
    log10_TPM = log10(TPM + 1e-6)
  )

# Prefer mixed model when site/depth structure is available.
fit <- lmer(log10_TPM ~ habitat + depth + (1 | site), data = dat)
anova_table <- anova(fit)
emm <- emmeans(fit, pairwise ~ habitat, adjust = "BH")

readr::write_csv(as.data.frame(anova_table), "habitat_mixed_model_anova.csv")
readr::write_csv(as.data.frame(emm$contrasts), "habitat_pairwise_BH.csv")

# If sample size is too small for the mixed model, Kruskal-Wallis can be used
# as exploratory analysis, but report that it does not control site nesting.
kw <- kruskal.test(log10_TPM ~ habitat, data = dat)
capture.output(kw, file = "habitat_kruskal_exploratory.txt")
EOF
  log "Wrote statistics template: ${rfile}"
}

#######################################
# 9. Step dispatcher
#######################################

usage() {
  cat <<EOF
Usage:
  bash $0 <step>

Steps:
  init                  Create folders and metadata templates
  check                 Validate samples.tsv and input FASTQ paths
  versions              Record software versions
  download              Download SRA and convert to FASTQ
  qc                    FastQC/MultiQC raw reads
  fastp                 fastp cleaning
  qc-clean              FastQC/MultiQC clean reads
  assembly              MEGAHIT assembly
  quast                 QUAST assembly QC
  metaphlan             MetaPhlAn profiling
  binning               MetaWRAP binning
  refinement            MetaWRAP bin refinement
  bin-stats             Combine bin stats
  collect-bins          Collect refined bins
  mag-qc                CheckM2 + GUNC + quality filtering
  drep                  dRep dereplication without --ignoreGenomeQuality
  gtdb                  GTDB-Tk classification and marker trees
  genes                 Prodigal gene prediction
  catalog               Unified non-redundant gene catalog
  salmon                Salmon mapping to unified catalog
  coverm                MAG abundance with CoverM
  kofam                 KOfamScan annotation
  k08356-search         PF00384/K08356 candidate discovery
  k08356-hmms           Build separate AioA and IdrA HMMs
  k08356-dual           Search candidates against AioA and IdrA HMMs
  k08356-tree           DMSOR-wide reference tree
  neighborhood          Create neighborhood evidence table template
  evidence              Create clade evidence table template
  stats-template        Write statistics R template
  all                   Run a recommended dry-run sequence

Important:
  Default is DRY_RUN=1.
  Use DRY_RUN=0 only after checking paths, conda env names and databases.
EOF
}

main() {
  local step="${1:-}"
  case "${step}" in
    ""|-h|--help|help) usage; return 0 ;;
  esac
  ensure_project_dirs
  case "${step}" in
    init) init_project ;;
    check) check_inputs ;;
    versions) record_versions ;;
    download) download_sra ;;
    qc) qc_raw ;;
    fastp) fastp_clean ;;
    qc-clean) qc_clean ;;
    assembly) assembly_megahit ;;
    quast) assembly_quast ;;
    metaphlan) metaphlan_profile ;;
    binning) metawrap_binning ;;
    refinement) metawrap_refinement ;;
    bin-stats) combine_bin_stats ;;
    collect-bins) collect_refined_bins ;;
    mag-qc) mag_quality_filter ;;
    drep) drep_all_samples ;;
    gtdb) gtdbtk_classify_and_tree ;;
    genes) predict_genes ;;
    catalog) build_gene_catalog ;;
    salmon) salmon_gene_abundance ;;
    coverm) mag_abundance_coverm ;;
    kofam) kofamscan_annotation ;;
    k08356-search) k08356_candidate_search ;;
    k08356-hmms) build_aioa_idra_hmms ;;
    k08356-dual) aioa_idra_dual_hmmsearch ;;
    k08356-tree) k08356_reference_tree ;;
    neighborhood) gene_neighborhood_template ;;
    evidence) clade_assignment_evidence_template ;;
    stats-template) write_statistics_template ;;
    all)
      init_project
      check_inputs
      record_versions
      qc_raw
      fastp_clean
      qc_clean
      assembly_megahit
      assembly_quast
      metaphlan_profile
      metawrap_binning
      metawrap_refinement
      combine_bin_stats
      collect_refined_bins
      mag_quality_filter
      drep_all_samples
      gtdbtk_classify_and_tree
      predict_genes
      build_gene_catalog
      salmon_gene_abundance
      mag_abundance_coverm
      kofamscan_annotation
      k08356_candidate_search
      build_aioa_idra_hmms
      aioa_idra_dual_hmmsearch
      k08356_reference_tree
      gene_neighborhood_template
      clade_assignment_evidence_template
      write_statistics_template
      ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
