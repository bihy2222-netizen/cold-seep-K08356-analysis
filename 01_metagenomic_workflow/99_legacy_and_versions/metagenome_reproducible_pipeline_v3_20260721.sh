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
# Version status:
#   Metagenome main workflow v1.2 + K08356 classification framework v1.0.
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
#   sample  habitat  site  core  depth_layer  depth_top_cm  depth_bottom_cm
#   depth_mid_cm  r1  r2  srr
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
RIESKE_PF00355_HMM="${RIESKE_PF00355_HMM:-${PROJECT_DIR}/references/PF00355.hmm}"
DIHAEM_PF03150_HMM="${DIHAEM_PF03150_HMM:-${PROJECT_DIR}/references/PF03150.hmm}"
AIOA_SEED_FAA="${AIOA_SEED_FAA:-${PROJECT_DIR}/references/aioA_seed.faa}"
IDRA_SEED_FAA="${IDRA_SEED_FAA:-${PROJECT_DIR}/references/idrA_seed.faa}"
DMSOR_REFERENCE_FAA="${DMSOR_REFERENCE_FAA:-${PROJECT_DIR}/references/dmsor_family_references.faa}"

LOG_DIR="${PROJECT_DIR}/logs"

ensure_project_dirs() {
  mkdir -p "${PROJECT_DIR}"/{metadata,logs,00_raw/sra,00_raw/fastq,01_qc/fastqc_raw,01_qc/fastqc_clean,01_qc/fastp,02_clean,03_assembly,03_assembly_qc,04_metaphlan,05_binning,06_genes/{raw/CDS,raw/AA,raw/GFF,CDS,AA,GFF,catalog},07_abundance/{gene,MAG},08_annotation/{kofam/tmp,k08356,hmm,hmm_calibration,neighborhood},09_MAGs/{all_refined_bins,filtered_MAGs,dRep_out,checkm2,gunc,gtdbtk},10_statistics,reports,references}
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
  sample<TAB>habitat<TAB>site<TAB>core<TAB>depth_layer<TAB>depth_top_cm<TAB>depth_bottom_cm<TAB>depth_mid_cm<TAB>r1<TAB>r2<TAB>srr

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
  cat > "${PROJECT_DIR}/metadata/samples.v2_template.tsv" <<'EOF'
sample	habitat	site	core	depth_layer	depth_top_cm	depth_bottom_cm	depth_mid_cm	r1	r2	srr
NS_0-6	NS	NS	NS_core1	surface	0	6	3	/path/to/NS_0-6_R1.fastq.gz	/path/to/NS_0-6_R2.fastq.gz	NA
ES_2_0-6	ES	ES_2	ES_2_core1	surface	0	6	3	/path/to/ES_2_0-6_R1.fastq.gz	/path/to/ES_2_0-6_R2.fastq.gz	NA
EOF
  if [[ ! -f "${SAMPLES_TSV}" ]]; then
    cp "${PROJECT_DIR}/metadata/samples.v2_template.tsv" "${SAMPLES_TSV}"
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
      required = "sample habitat site core depth_layer depth_top_cm depth_bottom_cm depth_mid_cm r1 r2 srr"
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
  for cmd in fastqc multiqc fastp fasterq-dump vdb-validate megahit quast.py metaphlan metawrap prodigal cd-hit-est cd-hit salmon checkm2 gunc dRep gtdb-tk trimal iqtree2 iqtree hmmbuild hmmpress hmmsearch mafft seqkit coverm gawk; do
    run "printf '${cmd}\t' >> '${out}'; (${cmd} --version 2>&1 | head -n 1 || true) >> '${out}'"
  done
}

check_k08356_references() {
  local failed=0 f
  for f in "${PFAM00384_HMM}" "${RIESKE_PF00355_HMM}" "${DIHAEM_PF03150_HMM}" \
           "${AIOA_SEED_FAA}" "${IDRA_SEED_FAA}" "${DMSOR_REFERENCE_FAA}" \
           "${PROJECT_DIR}/metadata/reference_protein_manifest.tsv"; do
    if [[ ! -s "${f}" ]]; then echo "Missing/empty K08356 reference: ${f}" >&2; failed=1; fi
  done
  [[ "${failed}" == 0 ]] || return 1
  for f in "${AIOA_SEED_FAA}" "${IDRA_SEED_FAA}" "${DMSOR_REFERENCE_FAA}"; do
    awk '/^>/{id=$1; sub(/^>/,"",id); if(seen[id]++) dup=1} END{exit dup?1:0}' "${f}" || {
      echo "Duplicate FASTA IDs in ${f}" >&2; return 1;
    }
  done
  log "K08356 reference files exist and FASTA identifiers are unique"
}

#######################################
# 3. Download and read QC
#######################################

download_sra() {
  activate_env sratools
  require_samples
  sample_rows | while IFS=$'\t' read -r sample habitat site core depth_layer depth_top_cm depth_bottom_cm depth_mid_cm r1 r2 srr; do
    [[ "${srr}" == "NA" || -z "${srr}" ]] && continue
    run "prefetch '${srr}' -O '${PROJECT_DIR}/00_raw/sra' > '${LOG_DIR}/${sample}.prefetch.log' 2>&1"
    run "vdb-validate '${PROJECT_DIR}/00_raw/sra/${srr}' > '${LOG_DIR}/${sample}.vdb_validate.log' 2>&1"
    run "fasterq-dump '${PROJECT_DIR}/00_raw/sra/${srr}' --split-files --threads 16 --outdir '${PROJECT_DIR}/00_raw/fastq' > '${LOG_DIR}/${sample}.fasterq_dump.log' 2>&1"
    run "gzip -f '${PROJECT_DIR}/00_raw/fastq/${srr}_1.fastq' '${PROJECT_DIR}/00_raw/fastq/${srr}_2.fastq'"
    run "ln -sf '${PROJECT_DIR}/00_raw/fastq/${srr}_1.fastq.gz' '${PROJECT_DIR}/00_raw/fastq/${sample}_R1.fastq.gz'"
    run "ln -sf '${PROJECT_DIR}/00_raw/fastq/${srr}_2.fastq.gz' '${PROJECT_DIR}/00_raw/fastq/${sample}_R2.fastq.gz'"
  done
  run "printf 'sample\\tr1\\tr2\\n' > '${PROJECT_DIR}/metadata/downloaded_fastq_paths.tsv'; awk 'BEGIN{FS=OFS=\"\\t\"} NR>1 && \$11 != \"NA\" {print \$1, \"'${PROJECT_DIR}'/00_raw/fastq/\"\$1\"_R1.fastq.gz\", \"'${PROJECT_DIR}'/00_raw/fastq/\"\$1\"_R2.fastq.gz\"}' '${SAMPLES_TSV}' >> '${PROJECT_DIR}/metadata/downloaded_fastq_paths.tsv'"
  run "md5sum '${PROJECT_DIR}'/00_raw/fastq/*.fastq* > '${PROJECT_DIR}/metadata/raw_fastq.md5' || true"
}

qc_raw() {
  activate_env fastqc
  require_samples
  sample_rows | while IFS=$'\t' read -r sample habitat site core depth_layer depth_top_cm depth_bottom_cm depth_mid_cm r1 r2 srr; do
    if [[ "${r1}" == "NA" && "${srr}" != "NA" ]]; then
      r1="${PROJECT_DIR}/00_raw/fastq/${sample}_R1.fastq.gz"
      r2="${PROJECT_DIR}/00_raw/fastq/${sample}_R2.fastq.gz"
    fi
    [[ "${r1}" == "NA" || "${r2}" == "NA" ]] && continue
    run "fastqc '${r1}' '${r2}' -o '${PROJECT_DIR}/01_qc/fastqc_raw' -t 2 > '${LOG_DIR}/${sample}.fastqc_raw.log' 2>&1"
  done
  activate_env multiqc
  run "multiqc '${PROJECT_DIR}/01_qc/fastqc_raw' -o '${PROJECT_DIR}/01_qc' > '${LOG_DIR}/multiqc_raw.log' 2>&1"
}

fastp_clean() {
  activate_env fastp
  require_samples
  sample_rows | while IFS=$'\t' read -r sample habitat site core depth_layer depth_top_cm depth_bottom_cm depth_mid_cm r1 r2 srr; do
    if [[ "${r1}" == "NA" && "${srr}" != "NA" ]]; then
      r1="${PROJECT_DIR}/00_raw/fastq/${sample}_R1.fastq.gz"
      r2="${PROJECT_DIR}/00_raw/fastq/${sample}_R2.fastq.gz"
    fi
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
  run "find '${PROJECT_DIR}/05_binning' -path '*/BIN_REFINEMENT/*bins.stats' -print > '${PROJECT_DIR}/09_MAGs/bin_refinement_stats_files.txt'"
  run "while read -r file; do [[ -f \"\$file\" ]] || continue; sample_dir=\$(basename \"\$(dirname \"\$(dirname \"\$file\")\")\"); sample_name=\${sample_dir%_bins}; awk -v sample=\"\$sample_name\" 'NR > 1 {gsub(/^bin\\./, \"\", \$1); printf \"%s\\t%s\", sample, \$1; for (i = 2; i <= NF; i++) printf \"\\t%s\", \$i; printf \"\\n\"}' OFS='\\t' \"\$file\" >> '${out}'; done < '${PROJECT_DIR}/09_MAGs/bin_refinement_stats_files.txt'"
}

collect_refined_bins() {
  # Dry-run this first. It copies bins and prefixes sample names to avoid ID collisions.
  sample_col sample | while read -r sample; do
    run "find '${PROJECT_DIR}/05_binning/${sample}_bins/BIN_REFINEMENT' -maxdepth 3 -type d -print > '${PROJECT_DIR}/09_MAGs/${sample}.BIN_REFINEMENT.directories.txt'"
    run "find '${PROJECT_DIR}/05_binning/${sample}_bins/BIN_REFINEMENT' -type f \\( -name '*.fa' -o -name '*.fasta' \\) ! -path '*/work_files/*' -print | while read -r bin; do base=\$(basename \"\$bin\"); base=\${base%.fa}; base=\${base%.fasta}; cp \"\$bin\" '${PROJECT_DIR}/09_MAGs/all_refined_bins/${sample}_'\"\$base\"'.fasta'; done"
  done
}

mag_quality_filter() {
  activate_env checkm2
  run "checkm2 predict --input '${PROJECT_DIR}/09_MAGs/all_refined_bins' --output-directory '${PROJECT_DIR}/09_MAGs/checkm2' --threads '${THREADS}' --extension fasta > '${LOG_DIR}/checkm2.log' 2>&1"

  activate_env gunc
  run "gunc run --input_dir '${PROJECT_DIR}/09_MAGs/all_refined_bins' --file_suffix .fasta --out_dir '${PROJECT_DIR}/09_MAGs/gunc' --threads 32 > '${LOG_DIR}/gunc.log' 2>&1"

  # Build a joint CheckM2 + GUNC table. GUNC failures are marked, not silently
  # deleted, because key K08356 MAGs may need manual chimera review.
  run "awk 'BEGIN{FS=OFS=\"\\t\"} FNR==NR {if(NR==1){for(i=1;i<=NF;i++) gh[\$i]=i; next} mag=\$1; gunc_pass[mag]=(\"pass\"); css[mag]=\"NA\"; rrs[mag]=\"NA\"; contam_portion[mag]=\"NA\"; for(k in gh){lk=tolower(k); if(lk ~ /pass/ || lk ~ /status/) gunc_pass[mag]=\$(gh[k]); if(lk==\"css\") css[mag]=\$(gh[k]); if(lk==\"rrs\") rrs[mag]=\$(gh[k]); if(lk ~ /contamination/) contam_portion[mag]=\$(gh[k])} next} FNR==1 {for(i=1;i<=NF;i++) h[\$i]=i; print \"MAG\", \"completeness\", \"contamination\", \"quality_score\", \"GUNC_pass\", \"CSS\", \"RRS\", \"contamination_portion\", \"final_keep\", \"exclusion_reason\"; next} FNR>1 {mag=\$1; comp=\$(h[\"Completeness\"]); con=\$(h[\"Contamination\"]); q=comp-5*con; gp=(mag in gunc_pass)?gunc_pass[mag]:\"not_available\"; reason=\"pass\"; keep=\"yes\"; if(comp < 50 || con > 10){keep=\"no\"; reason=\"low_quality\"} else if(gp !~ /pass|PASS|TRUE|true|not_available/){keep=\"review\"; reason=\"suspected_chimeric_MAG\"} print mag, comp, con, q, gp, ((mag in css)?css[mag]:\"NA\"), ((mag in rrs)?rrs[mag]:\"NA\"), ((mag in contam_portion)?contam_portion[mag]:\"NA\"), keep, reason}' '${PROJECT_DIR}/09_MAGs/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv' '${PROJECT_DIR}/09_MAGs/checkm2/quality_report.tsv' > '${PROJECT_DIR}/09_MAGs/MAG_quality_CheckM2_GUNC.tsv'"
  run "awk 'BEGIN{FS=OFS=\"\\t\"} NR>1 && (\$9 == \"yes\" || \$9 == \"review\") {print \$1}' '${PROJECT_DIR}/09_MAGs/MAG_quality_CheckM2_GUNC.tsv' | while read -r mag; do if [[ -f '${PROJECT_DIR}/09_MAGs/all_refined_bins/'\"\$mag\"'.fasta' ]]; then cp '${PROJECT_DIR}/09_MAGs/all_refined_bins/'\"\$mag\"'.fasta' '${PROJECT_DIR}/09_MAGs/filtered_MAGs/'; else echo \"Missing MAG fasta after QC: \$mag\" >&2; exit 1; fi; done"
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
    run "prodigal -i '${PROJECT_DIR}/03_assembly/${sample}/final.contigs.fa' -d '${PROJECT_DIR}/06_genes/raw/CDS/${sample}_cds.fa' -a '${PROJECT_DIR}/06_genes/raw/AA/${sample}_AA.faa' -f gff -o '${PROJECT_DIR}/06_genes/raw/GFF/${sample}.gff' -p meta -m > '${LOG_DIR}/${sample}.prodigal.log' 2>&1"
    # MEGAHIT/Prodigal identifiers recur between assemblies. Prefix sample IDs
    # before pooling so HMM hits, abundance and neighborhoods remain traceable.
    run "awk -v s='${sample}' '/^>/{sub(/^>/,\">\" s \"|\")}1' '${PROJECT_DIR}/06_genes/raw/CDS/${sample}_cds.fa' > '${PROJECT_DIR}/06_genes/CDS/${sample}_cds.fa'"
    run "awk -v s='${sample}' '/^>/{sub(/^>/,\">\" s \"|\")}1' '${PROJECT_DIR}/06_genes/raw/AA/${sample}_AA.faa' > '${PROJECT_DIR}/06_genes/AA/${sample}_AA.faa'"
    run "awk -v s='${sample}' 'BEGIN{FS=OFS=\"\\t\"} /^#/ {print; next} {\$1=s \"|\" \$1; n=split(\$9,a,\";\"); for(i=1;i<=n;i++){if(a[i] ~ /^(ID|Parent)=/){split(a[i],b,\"=\"); a[i]=b[1] \"=\" s \"|\" b[2]}} \$9=a[1]; for(i=2;i<=n;i++) \$9=\$9 \";\" a[i]; print}' '${PROJECT_DIR}/06_genes/raw/GFF/${sample}.gff' > '${PROJECT_DIR}/06_genes/GFF/${sample}.gff'"
  done
}

build_gene_catalog() {
  # Unified non-redundant catalog for cross-habitat comparison.
  # Critical design: cluster CDS only, then extract matching protein IDs from
  # Prodigal AA files. Do not independently cluster CDS and AA, because that
  # can choose different representatives and break TPM-to-protein annotation.
  activate_env cd-hit
  run "cat '${PROJECT_DIR}'/06_genes/CDS/*_cds.fa > '${PROJECT_DIR}/06_genes/catalog/all_samples_cds.fa'"
  run "cd-hit-est -i '${PROJECT_DIR}/06_genes/catalog/all_samples_cds.fa' -o '${PROJECT_DIR}/06_genes/catalog/all_samples_cds.nr95.fa' -c 0.95 -aS 0.9 -G 0 -M 0 -g 1 -T '${THREADS}' > '${LOG_DIR}/cdhit_catalog_cds.log' 2>&1"
  run "cat '${PROJECT_DIR}'/06_genes/AA/*_AA.faa > '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.faa'"

  activate_env seqkit
  run "seqkit seq -n '${PROJECT_DIR}/06_genes/catalog/all_samples_cds.nr95.fa' | sed 's/[[:space:]].*$//' > '${PROJECT_DIR}/06_genes/catalog/representative_cds_ids.txt'"
  run "seqkit grep -n -f '${PROJECT_DIR}/06_genes/catalog/representative_cds_ids.txt' '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.faa' > '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.from_cds_nr95.faa'"

  # Map each CD-HIT cluster representative to its member genes.
  run "awk 'BEGIN{OFS=\"\\t\"; print \"catalog_gene_id\", \"representative_CDS_id\", \"representative_protein_id\", \"member_gene_id\", \"sample\", \"contig\", \"start\", \"end\", \"strand\"} /^>Cluster/ {cluster=\$2; rep=\"\"; next} /^[0-9]/ {line=\$0; id=line; sub(/^.*>/, \"\", id); sub(/\\.\\.\\..*$/, \"\", id); if (line ~ /\\*/) rep=id; members[cluster]=members[cluster] id \"\\n\"; reps[cluster]=rep} END {for (c in members) {split(members[c], a, \"\\n\"); for (i in a) if (a[i] != \"\") {split(a[i],z,\"|\"); sample=z[1]; print \"catalog_\" c, reps[c], reps[c], a[i], sample, \"NA\", \"NA\", \"NA\", \"NA\"}}}' '${PROJECT_DIR}/06_genes/catalog/all_samples_cds.nr95.fa.clstr' > '${PROJECT_DIR}/06_genes/catalog/catalog_cds_protein_member_map.tsv'"

  # Add genomic coordinates from GFF when IDs match Prodigal IDs.
  run "cat '${PROJECT_DIR}'/06_genes/GFF/*.gff > '${PROJECT_DIR}/06_genes/catalog/all_samples.unique_ids.gff'"
  run "awk 'BEGIN{FS=OFS=\"\\t\"} FNR==NR && \$0 !~ /^#/ {id=\"\"; split(\$9,a,\";\"); for(i in a){if(a[i] ~ /^ID=/){id=a[i]; sub(/^ID=/,\"\",id)}} if(id != \"\") coord[id]=\$1 OFS \$4 OFS \$5 OFS \$7; next} FNR==1 {print; next} {if(\$4 in coord){split(coord[\$4],c,OFS); \$6=c[1]; \$7=c[2]; \$8=c[3]; \$9=c[4]} print}' '${PROJECT_DIR}/06_genes/catalog/all_samples.unique_ids.gff' '${PROJECT_DIR}/06_genes/catalog/catalog_cds_protein_member_map.tsv' > '${PROJECT_DIR}/06_genes/catalog/catalog_cds_protein_member_map.with_coords.tsv'"
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
  run "'${KOFAM_EXEC}' --cpu '${THREADS}' -f detail-tsv --profile '${KOFAM_PROFILES}' --ko-list '${KOFAM_KO_LIST}' --tmp-dir '${PROJECT_DIR}/08_annotation/kofam/tmp' -o '${PROJECT_DIR}/08_annotation/kofam/all_samples_AA.from_cds_nr95.kofamscan.tsv' '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.from_cds_nr95.faa' > '${LOG_DIR}/kofamscan.log' 2>&1"
}

k08356_candidate_search() {
  activate_env hmmer
  run "hmmpress '${PFAM00384_HMM}' || true"
  run "hmmsearch --cpu '${THREADS}' --domtblout '${PROJECT_DIR}/08_annotation/k08356/PF00384.domtblout' --tblout '${PROJECT_DIR}/08_annotation/k08356/PF00384.tblout' '${PFAM00384_HMM}' '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.from_cds_nr95.faa' > '${LOG_DIR}/PF00384.hmmsearch.log' 2>&1"

  # KOfam K08356 candidates, if KOfam output includes KO identifiers.
  run "awk 'BEGIN{FS=OFS=\"\\t\"; print \"candidate_id\",\"source\",\"significant\",\"raw_line\"} /K08356/ {sig=(\$1==\"*\"?\"yes\":\"no\"); id=(\$1==\"*\"?\$2:\$1); if(id!=\"\") print id,\"KOfam_K08356\",sig,\$0}' '${PROJECT_DIR}/08_annotation/kofam/all_samples_AA.from_cds_nr95.kofamscan.tsv' > '${PROJECT_DIR}/08_annotation/k08356/K08356.filtered_hits.tsv'"

  # Parse PF00384 HMMER domain hits. Coverage is domain alignment length / target length.
  run "awk 'BEGIN{FS=OFS=\"\\t\"; print \"candidate_id\", \"source\", \"full_evalue\", \"full_score\", \"domain_ievalue\", \"domain_score\", \"target_len\", \"domain_cov\"} \$0 !~ /^#/ {cov=(\$19-\$18+1)/\$3; if(\$7 <= 1e-5 && cov >= 0.30) print \$1, \"PF00384\", \$7, \$8, \$13, \$14, \$3, cov}' '${PROJECT_DIR}/08_annotation/k08356/PF00384.domtblout' > '${PROJECT_DIR}/08_annotation/k08356/PF00384.filtered_hits.tsv'"

  run "awk 'BEGIN{FS=OFS=\"\\t\"} FNR==1 {next} {print \$1}' '${PROJECT_DIR}/08_annotation/k08356/PF00384.filtered_hits.tsv' '${PROJECT_DIR}/08_annotation/k08356/K08356.filtered_hits.tsv' | sort -u > '${PROJECT_DIR}/08_annotation/k08356/k08356_candidate_ids.txt'"
  activate_env seqkit
  run "seqkit grep -n -f '${PROJECT_DIR}/08_annotation/k08356/k08356_candidate_ids.txt' '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.from_cds_nr95.faa' > '${PROJECT_DIR}/08_annotation/k08356/k08356_candidates.faa'"
  run "awk 'BEGIN{FS=OFS=\"\\t\"; print \"candidate_id\", \"source\", \"note\"} {print \$1, \"PF00384_or_K08356\", \"Extracted to k08356_candidates.faa\"}' '${PROJECT_DIR}/08_annotation/k08356/k08356_candidate_ids.txt' > '${PROJECT_DIR}/08_annotation/k08356/k08356_candidates.tsv'"
  run "seqkit stats '${PROJECT_DIR}/08_annotation/k08356/k08356_candidates.faa' > '${PROJECT_DIR}/08_annotation/k08356/k08356_candidates.seqkit_stats.tsv'"
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
  # Score independent curated references as the calibration set. Training-only
  # sequences must be marked separately in reference_protein_manifest.tsv.
  run "hmmsearch --cpu '${THREADS}' --domtblout '${PROJECT_DIR}/08_annotation/hmm_calibration/idrA.references.domtblout' '${PROJECT_DIR}/08_annotation/hmm/idrA.hmm' '${DMSOR_REFERENCE_FAA}' > '${LOG_DIR}/idrA.references.hmmsearch.log' 2>&1"
  run "hmmsearch --cpu '${THREADS}' --domtblout '${PROJECT_DIR}/08_annotation/hmm_calibration/aioA.references.domtblout' '${PROJECT_DIR}/08_annotation/hmm/aioA.hmm' '${DMSOR_REFERENCE_FAA}' > '${LOG_DIR}/aioA.references.hmmsearch.log' 2>&1"
}

k08356_calibrate_hmm_thresholds() {
  # Formal calibration should use independent positive/negative references and
  # cross-validation. This step creates an explicit threshold file instead of
  # hiding cutoffs inside the code. Edit values after inspecting reference scores.
  local out="${PROJECT_DIR}/08_annotation/hmm_calibration/aioA_idrA_thresholds.tsv"
  if [[ ! -f "${out}" ]]; then
    run "printf 'threshold_version\\tstatus\\tmin_domain_score\\tmin_target_cov\\tmin_model_cov\\tmin_score_delta\\tgray_zone_delta\\tcalibration_note\\n' > '${out}'"
    run "printf 'UNVALIDATED_v0.1\\tPLACEHOLDER\\tNA\\t0.50\\t0.50\\tNA\\tNA\\tDO_NOT_CLASSIFY: calibrate with independent AioA/IdrA positives, DMSOR negatives, tree placement and complete neighborhoods\\n' >> '${out}'"
  fi

  # Retain only the best-scoring domain per protein. Report both target and
  # model coverage; a short fragment can otherwise appear deceptively complete.
  for family in aioA idrA; do
    local label="AioA"; [[ "${family}" == "idrA" ]] && label="IdrA"
    for setname in candidate reference; do
      local src="${PROJECT_DIR}/08_annotation/k08356/${family}.domtblout"
      [[ "${setname}" == "reference" ]] && src="${PROJECT_DIR}/08_annotation/hmm_calibration/${family}.references.domtblout"
      local dst="${PROJECT_DIR}/08_annotation/hmm_calibration/${family}_${setname}_scores.tsv"
      run "awk -v fam='${label}' 'BEGIN{FS=OFS=\"\\t\"} \$0 !~ /^#/ {id=\$1; if(!(id in best) || \$14>best[id]){best[id]=\$14; line[id]=\$0}} END{print \"family\",\"sequence_id\",\"full_evalue\",\"full_score\",\"domain_ievalue\",\"domain_score\",\"target_len\",\"model_len\",\"hmm_from\",\"hmm_to\",\"ali_from\",\"ali_to\",\"target_cov\",\"model_cov\"; for(id in line){split(line[id],x,/ +/); tc=(x[19]-x[18]+1)/x[3]; mc=(x[17]-x[16]+1)/x[6]; print fam,x[1],x[7],x[8],x[13],x[14],x[3],x[6],x[16],x[17],x[18],x[19],tc,mc}}' '${src}' > '${dst}'"
    done
  done

  run "printf 'IMPORTANT: replace PLACEHOLDER only after joining *_reference_scores.tsv to metadata/reference_protein_manifest.tsv and validating against phylogeny plus idrABP1P2/aioBA neighborhoods. The published score 640 is not transferable to a newly built HMM.\\n' > '${PROJECT_DIR}/08_annotation/hmm_calibration/READ_BEFORE_CLASSIFICATION.txt'"
}

k08356_classify_by_hmm() {
  local threshold_file="${PROJECT_DIR}/08_annotation/hmm_calibration/aioA_idrA_thresholds.tsv"
  local out="${PROJECT_DIR}/08_annotation/k08356/k08356_dual_hmm_classification.tsv"
  if awk 'BEGIN{FS=\"\\t\"} NR==2 && (\$2==\"PLACEHOLDER\" || \$3==\"NA\") {exit 0} NR==2 {exit 1}' "${threshold_file}"; then
    echo "Refusing classification: HMM thresholds are still PLACEHOLDER/UNVALIDATED." >&2
    return 2
  fi
  run "gawk 'BEGIN{FS=OFS=\"\\t\"} FNR==NR && NR>1 {version=\$1; min_score=\$3; min_tcov=\$4; min_mcov=\$5; min_delta=\$6; gray=\$7; next} FILENAME==ARGV[2] && FNR>1 {as[\$2]=\$6; af[\$2]=\$4; at[\$2]=\$13; am[\$2]=\$14; ids[\$2]=1; next} FILENAME==ARGV[3] && FNR>1 {is[\$2]=\$6; inf[\$2]=\$4; it[\$2]=\$13; im[\$2]=\$14; ids[\$2]=1; next} END{print \"sequence_id\",\"AioA_full_score\",\"AioA_domain_score\",\"AioA_target_cov\",\"AioA_model_cov\",\"IdrA_full_score\",\"IdrA_domain_score\",\"IdrA_target_cov\",\"IdrA_model_cov\",\"score_delta\",\"HMM_class\",\"threshold_version\"; for(id in ids){a=as[id]+0; i=is[id]+0; d=a-i; cls=\"unknown_by_HMM\"; if(a>=min_score && at[id]>=min_tcov && am[id]>=min_mcov && d>=min_delta) cls=\"AioA_preliminary\"; else if(i>=min_score && it[id]>=min_tcov && im[id]>=min_mcov && -d>=min_delta) cls=\"IdrA_preliminary\"; else if((a>=min_score || i>=min_score) && d<gray && d>-gray) cls=\"gray_zone_unknown\"; print id,af[id]+0,a,at[id]+0,am[id]+0,inf[id]+0,i,it[id]+0,im[id]+0,d,cls,version}}' '${threshold_file}' '${PROJECT_DIR}/08_annotation/hmm_calibration/aioA_candidate_scores.tsv' '${PROJECT_DIR}/08_annotation/hmm_calibration/idrA_candidate_scores.tsv' > '${out}'"
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
  # Extract +/-10 ORF neighborhoods from Prodigal GFF for K08356 candidates.
  # This creates a real neighborhood scaffold; functional labels for neighbor
  # ORFs should be filled by joining KOfam/DRAM/Prokka annotations afterward.
  local out="${PROJECT_DIR}/08_annotation/neighborhood/k08356_gene_neighborhood.tsv"
  run "awk 'BEGIN{FS=OFS=\"\\t\"} FNR==NR {targets[\$1]=1; next} \$0 !~ /^#/ && \$3 == \"CDS\" {id=\"\"; split(\$9,a,\";\"); for(i in a){if(a[i] ~ /^ID=/){id=a[i]; sub(/^ID=/,\"\",id)}} if(id==\"\") next; key=FILENAME \"|\" \$1; n[key]++; idx=key SUBSEP n[key]; gene[idx]=id; start[idx]=\$4; end[idx]=\$5; strand[idx]=\$7; contig[idx]=\$1; file_idx[id]=key; pos_idx[id]=n[key]} END{print \"target_id\", \"contig\", \"contig_length\", \"target_start\", \"target_end\", \"distance_to_left_edge\", \"distance_to_right_edge\", \"upstream_10_orfs\", \"downstream_10_orfs\", \"aioB_or_idrB\", \"P1_P2\", \"electron_transfer_gene\", \"arsenic_resistance_gene\", \"transposase_integrase\", \"neighborhood_class\"; for(t in targets){key=file_idx[t]; p=pos_idx[t]; if(key==\"\" || p==\"\"){print t,\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"NA\",\"not_found_in_GFF\"; continue} up=\"\"; down=\"\"; for(i=p-10;i<p;i++){if(i>=1){idx=key SUBSEP i; up=(up==\"\"?gene[idx]:up\",\"gene[idx])}} for(i=p+1;i<=p+10;i++){idx=key SUBSEP i; if(gene[idx] != \"\") down=(down==\"\"?gene[idx]:down\",\"gene[idx])} idx=key SUBSEP p; edge=\"requires_contig_length\"; print t, contig[idx], \"NA\", start[idx], end[idx], \"NA\", \"NA\", up, down, \"to_annotate\", \"to_annotate\", \"to_annotate\", \"to_annotate\", \"to_annotate\", edge}}' '${PROJECT_DIR}/08_annotation/k08356/k08356_candidate_ids.txt' '${PROJECT_DIR}'/06_genes/GFF/*.gff > '${out}'"
}

annotate_neighborhood_markers() {
  # Search the complete, uniquely named protein collection for the two marker
  # families used in the ISMEJ gene-neighborhood comparison. These hit lists
  # are joined to the +/-10 ORF table during manual/ scripted evidence review.
  activate_env hmmer
  run "hmmsearch --cpu '${THREADS}' --domtblout '${PROJECT_DIR}/08_annotation/neighborhood/PF00355_Rieske.domtblout' '${RIESKE_PF00355_HMM}' '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.faa' > '${LOG_DIR}/PF00355_Rieske.log' 2>&1"
  run "hmmsearch --cpu '${THREADS}' --domtblout '${PROJECT_DIR}/08_annotation/neighborhood/PF03150_dihaem.domtblout' '${DIHAEM_PF03150_HMM}' '${PROJECT_DIR}/06_genes/catalog/all_samples_AA.faa' > '${LOG_DIR}/PF03150_dihaem.log' 2>&1"
  run "awk '\$0!~/^#/ && \$13<=1e-5 {print \$1}' '${PROJECT_DIR}/08_annotation/neighborhood/PF00355_Rieske.domtblout' | sort -u > '${PROJECT_DIR}/08_annotation/neighborhood/PF00355_Rieske.ids'"
  run "awk '\$0!~/^#/ && \$13<=1e-5 {print \$1}' '${PROJECT_DIR}/08_annotation/neighborhood/PF03150_dihaem.domtblout' | sort -u > '${PROJECT_DIR}/08_annotation/neighborhood/PF03150_dihaem.ids'"
  run "printf 'Interpretation gate:\n- Rieske plus >=2 adjacent PF03150 hits supports an IdrA-like idrABP1P2 neighborhood.\n- Rieske without adjacent PF03150 supports an AioA-like aioBA neighborhood.\n- Always review contig boundaries; absence on a truncated contig is not biological absence.\n' > '${PROJECT_DIR}/08_annotation/neighborhood/NEIGHBORHOOD_INTERPRETATION.txt'"
}

clade_assignment_evidence_template() {
  # Final class must combine HMM, tree and neighborhood evidence.
  # Classes: canonical_AioA_associated,
  # AioA_proximal_uncharacterized_DMSOR, IdrA_associated,
  # IdrA_proximal_uncharacterized_DMSOR. HMM alone never assigns function.
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

packages <- c("tidyverse", "lme4", "lmerTest", "emmeans", "effectsize", "performance")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)
library(performance)

abund_file <- "../07_abundance/gene/k08356_clade_abundance.tsv"
meta_file <- "../metadata/samples.tsv"

abund <- readr::read_tsv(abund_file, show_col_types = FALSE)
meta <- readr::read_tsv(meta_file, show_col_types = FALSE)

dat <- abund %>%
  left_join(meta, by = "sample") %>%
  mutate(
    habitat = factor(habitat),
    site = factor(site),
    core = factor(core),
    depth_mid_cm = as.numeric(depth_mid_cm),
    log10_TPM = log10(TPM + 1)
  )

# Prefer mixed model when site/depth structure is available.
fit <- lmer(log10_TPM ~ habitat + depth_mid_cm + (1 | site), data = dat)
anova_table <- anova(fit)
emm <- emmeans(fit, pairwise ~ habitat, adjust = "BH")

group_n <- dat %>% count(habitat, name = "n")
readr::write_csv(as.data.frame(anova_table), "habitat_mixed_model_anova.csv")
readr::write_csv(as.data.frame(confint(emm$contrasts)), "habitat_pairwise_BH_95CI.csv")
readr::write_csv(group_n, "habitat_sample_size.csv")
capture.output(performance::check_model(fit), file = "model_diagnostics.txt")
capture.output(performance::r2_nakagawa(fit), file = "mixed_model_R2.txt")
capture.output(effectsize::eta_squared(anova_table), file = "effect_size_eta_squared.txt")
capture.output(isSingular(fit), file = "singular_fit_check.txt")

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
  k08356-check-refs     Validate HMM/reference inputs and unique FASTA IDs
  k08356-search         PF00384/K08356 candidate discovery
  k08356-hmms           Build separate AioA and IdrA HMMs
  k08356-dual           Search candidates against AioA and IdrA HMMs
  k08356-calibrate      Write/parse HMM threshold calibration tables
  k08356-classify       Classify candidates by dual-HMM score table
  k08356-tree           DMSOR-wide reference tree
  neighborhood          Create neighborhood evidence table template
  neighborhood-markers  Search PF00355 Rieske and PF03150 di-haem markers
  evidence              Create clade evidence table template
  stats-template        Write statistics R template
  dry-run-plan          Print a recommended dry-run sequence
  all                   Alias of dry-run-plan; forbidden when DRY_RUN=0

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
    k08356-check-refs) check_k08356_references ;;
    k08356-search) k08356_candidate_search ;;
    k08356-hmms) build_aioa_idra_hmms ;;
    k08356-dual) aioa_idra_dual_hmmsearch ;;
    k08356-calibrate) k08356_calibrate_hmm_thresholds ;;
    k08356-classify) k08356_classify_by_hmm ;;
    k08356-tree) k08356_reference_tree ;;
    neighborhood) gene_neighborhood_template ;;
    neighborhood-markers) annotate_neighborhood_markers ;;
    evidence) clade_assignment_evidence_template ;;
    stats-template) write_statistics_template ;;
    dry-run-plan|all)
      if [[ "${DRY_RUN}" != "1" ]]; then
        cat >&2 <<'EOF'
Refusing to run the entire workflow with DRY_RUN=0.

This pipeline intentionally requires stepwise execution and inspection:
  check -> qc -> fastp -> qc-clean -> assembly -> quast -> ...

Use DRY_RUN=1 for all/dry-run-plan, then run individual steps with DRY_RUN=0.
EOF
        exit 1
      fi
      init_project
      cat <<EOF
Recommended stepwise execution after editing samples.tsv and reference paths:
  bash $0 check
  bash $0 versions
  bash $0 qc
  bash $0 fastp
  bash $0 qc-clean
  bash $0 assembly
  bash $0 quast
  bash $0 metaphlan
  bash $0 binning
  bash $0 refinement
  bash $0 bin-stats
  bash $0 collect-bins
  bash $0 mag-qc
  bash $0 drep
  bash $0 gtdb
  bash $0 genes
  bash $0 catalog
  bash $0 salmon
  bash $0 coverm
	  bash $0 kofam
	  bash $0 k08356-check-refs
  bash $0 k08356-search
  bash $0 k08356-hmms
  bash $0 k08356-dual
  bash $0 k08356-calibrate
  bash $0 k08356-classify
  bash $0 k08356-tree
	  bash $0 neighborhood
	  bash $0 neighborhood-markers
  bash $0 evidence
  bash $0 stats-template
EOF
      ;;
    *) usage; exit 1 ;;
  esac
}

main "$@"
