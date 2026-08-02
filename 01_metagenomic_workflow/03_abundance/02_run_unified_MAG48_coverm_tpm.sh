#!/usr/bin/env bash
set -euo pipefail

# Build one exact 48-MAG reference and profile all 56 samples against it with
# identical CoverM settings. RUN_COVERM=0 prepares and validates inputs only.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ROOT_DIR="${ROOT_DIR:-/home/ps/ps1/data/bihongyu/cold_seep/illu/binning}"
RUN_DIR="${RUN_DIR:-${ROOT_DIR}/K08356_four_group_annotation/unified_MAG48_coverm_tpm_20260802}"
MAG_MANIFEST="${MAG_MANIFEST:-${REPO_ROOT}/01_metagenomic_workflow/04_metabolism/METABOLIC_group_derep49/manifests/MAG49_input_manifest.tsv}"
SAMPLE_HABITAT_CSV="${SAMPLE_HABITAT_CSV:-${REPO_ROOT}/02_figure_scripts/10_k08356_metabolism/sample_group_corrected.csv}"
OLD_DEREP_DIR="${OLD_DEREP_DIR:-${ROOT_DIR}/all_50_10_dRep/dereplicated_genomes}"
COVERM_BIN="${COVERM_BIN:-/home/ps/anaconda3/envs/coverm/bin/coverm}"
THREADS="${THREADS:-100}"
RUN_COVERM="${RUN_COVERM:-0}"

REFERENCE_DIR="${RUN_DIR}/reference_MAG48"
TPM_DIR="${RUN_DIR}/tpm_per_sample"
MANIFEST_DIR="${RUN_DIR}/manifests"
LOG_DIR="${RUN_DIR}/logs"
SAMPLE_TABLE="${MANIFEST_DIR}/sample56_reads_habitat.tsv"
REFERENCE_MANIFEST="${MANIFEST_DIR}/MAG48_reference_manifest.tsv"
SETTINGS_MANIFEST="${MANIFEST_DIR}/coverm_settings.tsv"

mkdir -p "$REFERENCE_DIR" "$TPM_DIR" "$MANIFEST_DIR" "$LOG_DIR"

for path in "$MAG_MANIFEST" "$SAMPLE_HABITAT_CSV" \
  "${OLD_DEREP_DIR}/coverm1.sh" "${OLD_DEREP_DIR}/coverm2.sh"; do
  [[ -s "$path" ]] || { echo "Error: missing input: $path" >&2; exit 1; }
done
[[ -x "$COVERM_BIN" ]] || {
  echo "Error: CoverM is not executable: $COVERM_BIN" >&2
  exit 1
}

tmp_dir="$(mktemp -d "${RUN_DIR}/prepare.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT

awk -F ',' '
  BEGIN { OFS = "\t" }
  NR > 1 {
    gsub(/"/, "", $1)
    gsub(/"/, "", $2)
    if ($1 != "" && $2 != "") print $1, $2
  }
' "$SAMPLE_HABITAT_CSV" > "${tmp_dir}/sample_habitat.tsv"

awk '
  BEGIN { OFS = "\t" }
  {
    r1 = ""; r2 = ""; out = ""
    for (i = 1; i <= NF; i++) {
      if ($i == "-1") r1 = $(i + 1)
      if ($i == "-2") r2 = $(i + 1)
      if ($i == "-o") out = $(i + 1)
    }
    if (r1 != "" && r2 != "" && out != "") {
      sample = out
      sub(/^bin_tpm\//, "", sample)
      sub(/_tpm\.txt$/, "", sample)
      print sample, r1, r2
    }
  }
' "${OLD_DEREP_DIR}/coverm1.sh" "${OLD_DEREP_DIR}/coverm2.sh" |
  sort -u > "${tmp_dir}/sample_reads.tsv"

{
  printf 'sample\thabitat\tr1\tr2\n'
  awk -F '\t' '
    BEGIN { OFS = "\t" }
    NR == FNR { habitat[$1] = $2; next }
    {
      if (!($1 in habitat)) {
        print "Missing habitat for sample " $1 > "/dev/stderr"
        exit 2
      }
      print $1, habitat[$1], $2, $3
    }
  ' "${tmp_dir}/sample_habitat.tsv" "${tmp_dir}/sample_reads.tsv" | sort
} > "$SAMPLE_TABLE"

sample_count="$(awk 'NR > 1 {print $1}' "$SAMPLE_TABLE" | sort -u | wc -l | tr -d ' ')"
sample_rows="$(awk 'NR > 1 {n++} END {print n+0}' "$SAMPLE_TABLE")"
[[ "$sample_count" == 56 && "$sample_rows" == 56 ]] || {
  echo "Error: expected 56 unique samples, found ${sample_count}/${sample_rows}." >&2
  exit 1
}

while IFS=$'\t' read -r sample habitat r1 r2; do
  [[ "$sample" == "sample" ]] && continue
  [[ -s "$r1" ]] || { echo "Error: missing R1 for $sample: $r1" >&2; exit 1; }
  [[ -s "$r2" ]] || { echo "Error: missing R2 for $sample: $r2" >&2; exit 1; }
done < "$SAMPLE_TABLE"

awk -F '\t' '
  BEGIN { OFS = "\t" }
  NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    print "source_habitat", "MAG", "source_fasta"
    next
  }
  {
    source = $(col["source_faa"])
    sub(/\.cds\.faa$/, ".fasta", source)
    print $(col["group"]), $(col["bin_id"]), source
  }
' "$MAG_MANIFEST" > "${tmp_dir}/MAG48_source.tsv"

{
  printf 'source_habitat\tMAG\tsource_fasta\tsha256\treference_symlink\n'
  while IFS=$'\t' read -r habitat MAG source_fasta; do
    [[ "$MAG" == "MAG" ]] && continue
    [[ -s "$source_fasta" ]] || {
      echo "Error: missing MAG genome FASTA: $source_fasta" >&2
      exit 1
    }
    link_path="${REFERENCE_DIR}/${MAG}.fasta"
    ln -sfn "$source_fasta" "$link_path"
    checksum="$(sha256sum "$source_fasta" | awk '{print $1}')"
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$habitat" "$MAG" "$source_fasta" "$checksum" "$link_path"
  done < "${tmp_dir}/MAG48_source.tsv"
} > "$REFERENCE_MANIFEST"

MAG_count="$(awk 'NR > 1 {print $2}' "$REFERENCE_MANIFEST" | sort -u | wc -l | tr -d ' ')"
MAG_rows="$(awk 'NR > 1 {n++} END {print n+0}' "$REFERENCE_MANIFEST")"
reference_count="$(find "$REFERENCE_DIR" -maxdepth 1 -type l -name '*.fasta' | wc -l | tr -d ' ')"
[[ "$MAG_count" == 48 && "$MAG_rows" == 48 && "$reference_count" == 48 ]] || {
  echo "Error: expected an exact 48-MAG reference; found ${MAG_count}/${MAG_rows}/${reference_count}." >&2
  exit 1
}

{
  printf 'parameter\tvalue\n'
  printf 'reference_MAG_count\t48\n'
  printf 'sample_count\t56\n'
  printf 'method\ttpm\n'
  printf 'min_covered_fraction\t0.1\n'
  printf 'min_read_percent_identity\t95\n'
  printf 'min_read_aligned_percent\t75\n'
  printf 'trim_min\t0.1\n'
  printf 'trim_max\t0.9\n'
  printf 'threads\t%s\n' "$THREADS"
  printf 'reference_directory\t%s\n' "$REFERENCE_DIR"
} > "$SETTINGS_MANIFEST"

printf '[%s] Validated common reference: 48 MAGs\n' "$(date '+%F %T')"
printf '[%s] Validated sample table: 56 samples\n' "$(date '+%F %T')"

if [[ "$RUN_COVERM" != 1 ]]; then
  printf '[%s] Preparation only. Set RUN_COVERM=1 to start profiling.\n' \
    "$(date '+%F %T')"
  exit 0
fi

while IFS=$'\t' read -r sample habitat r1 r2; do
  [[ "$sample" == "sample" ]] && continue
  out_file="${TPM_DIR}/${sample}_MAG48_tpm.tsv"
  log_file="${LOG_DIR}/${sample}.coverm.log"
  if [[ -s "$out_file" ]]; then
    printf '[%s] SKIP %s existing %s\n' "$(date '+%F %T')" "$sample" "$out_file"
    continue
  fi
  printf '[%s] RUN %s (%s) against common MAG48 reference\n' \
    "$(date '+%F %T')" "$sample" "$habitat"
  "$COVERM_BIN" genome \
    --genome-fasta-directory "$REFERENCE_DIR" \
    --genome-fasta-extension fasta \
    -1 "$r1" \
    -2 "$r2" \
    -m tpm \
    --min-covered-fraction 0.1 \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    --trim-min 0.1 \
    --trim-max 0.9 \
    -t "$THREADS" \
    -o "$out_file" > "$log_file" 2>&1
done < "$SAMPLE_TABLE"

RUN_DIR="$RUN_DIR" Rscript "${SCRIPT_DIR}/03_audit_unified_MAG48_tpm.R"
printf '[%s] Unified 56-sample x 48-MAG TPM workflow completed.\n' \
  "$(date '+%F %T')"
