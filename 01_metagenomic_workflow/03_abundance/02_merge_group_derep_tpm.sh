#!/usr/bin/env bash
set -euo pipefail

# Merge per-sample CoverM TPM files from the IS/AS/ES/NS group-derep MAG sets.
# The output is a tidy long table with explicit group and sample columns.

ROOT_DIR="${ROOT_DIR:-/home/ps/ps1/data/bihongyu/cold_seep/illu/binning}"
OUTPUT_DIR="${OUTPUT_DIR:-${ROOT_DIR}/group_tpm_merged}"
SAMPLE_TABLE="${SAMPLE_TABLE:-${ROOT_DIR}/group_coverm_samples.tsv}"
HABITAT_GROUPS=(IS AS ES NS)

group_genome_dir() {
  local group="$1"
  case "$group" in
    IS) printf '%s/IS/IS_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    AS) printf '%s/AS/AS_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    ES) printf '%s/ES/ES_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    NS) printf '%s/NS/NS_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    *) echo "Error: unknown group: $group" >&2; return 1 ;;
  esac
}

mkdir -p "$OUTPUT_DIR"
tmp_dir="$(mktemp -d "${OUTPUT_DIR}/.merge_group_tpm.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT

all_tmp="${tmp_dir}/all_groups_MAG_TPM_long.tsv"
summary_tmp="${tmp_dir}/all_groups_MAG_TPM_summary.tsv"
printf 'Group\tSample\tGenome\tTPM\tSource_file\n' > "$all_tmp"
printf 'Group\tTPM_files\tGenome_rows_per_file\tExpected_records\tObserved_records\tNonzero_TPM_rows\tTotal_TPM\tInput_directory\n' > "$summary_tmp"

expected_samples=""
if [[ -s "$SAMPLE_TABLE" ]]; then
  expected_samples="$(awk 'NF {n++} END {print n+0}' "$SAMPLE_TABLE")"
fi

for group in "${HABITAT_GROUPS[@]}"; do
  genome_dir="$(group_genome_dir "$group")"
  tpm_dir="${genome_dir}/bin_tpm"
  [[ -d "$tpm_dir" ]] || { echo "Error: TPM directory not found: $tpm_dir" >&2; exit 1; }

  mapfile -t tpm_files < <(find "$tpm_dir" -maxdepth 1 -type f -name '*_tpm.txt' -size +0c | sort)
  file_count="${#tpm_files[@]}"
  [[ "$file_count" -gt 0 ]] || { echo "Error: no non-empty TPM files in $tpm_dir" >&2; exit 1; }

  if [[ -n "$expected_samples" && "$file_count" -ne "$expected_samples" ]]; then
    echo "Error: group $group has $file_count TPM files; expected $expected_samples" >&2
    exit 1
  fi

  genome_count="$(find "$genome_dir" -maxdepth 1 -type f -name '*.fasta' | wc -l | tr -d ' ')"
  [[ "$genome_count" -gt 0 ]] || { echo "Error: no MAG FASTA files in $genome_dir" >&2; exit 1; }

  group_tmp="${tmp_dir}/${group}_MAG_TPM_long.tsv"
  printf 'Group\tSample\tGenome\tTPM\tSource_file\n' > "$group_tmp"

  for tpm_file in "${tpm_files[@]}"; do
    sample="$(basename "$tpm_file" _tpm.txt)"
    source_file="$(basename "$tpm_file")"
    record_count="$(awk 'END {print (NR > 0 ? NR - 1 : 0)}' "$tpm_file")"
    if [[ "$record_count" -ne "$genome_count" ]]; then
      echo "Error: $tpm_file has $record_count MAG rows; expected $genome_count" >&2
      exit 1
    fi

    awk -F '\t' -v OFS='\t' -v group="$group" -v sample="$sample" -v source="$source_file" '
      NR == 1 {
        if (NF != 2 || $1 != "Genome") {
          print "Error: unexpected CoverM header in " FILENAME > "/dev/stderr"
          exit 2
        }
        next
      }
      NF != 2 || $1 == "" || $2 == "" {
        print "Error: malformed TPM row " NR " in " FILENAME > "/dev/stderr"
        exit 3
      }
      {
        print group, sample, $1, $2, source
      }
    ' "$tpm_file" >> "$group_tmp"
  done

  tail -n +2 "$group_tmp" >> "$all_tmp"
  observed_records="$(awk 'END {print (NR > 0 ? NR - 1 : 0)}' "$group_tmp")"
  expected_records="$((file_count * genome_count))"
  [[ "$observed_records" -eq "$expected_records" ]] || {
    echo "Error: group $group merged $observed_records records; expected $expected_records" >&2
    exit 1
  }

  read -r nonzero_rows total_tpm < <(
    awk -F '\t' 'NR > 1 {if (($4 + 0) != 0) nonzero++; total += ($4 + 0)} END {printf "%d %.10g\n", nonzero+0, total+0}' "$group_tmp"
  )
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$group" "$file_count" "$genome_count" "$expected_records" \
    "$observed_records" "$nonzero_rows" "$total_tpm" "$tpm_dir" >> "$summary_tmp"

  mv "$group_tmp" "${OUTPUT_DIR}/${group}_MAG_TPM_long.tsv"
  printf '[%s] %s: %s files x %s MAGs = %s records\n' \
    "$(date '+%F %T')" "$group" "$file_count" "$genome_count" "$observed_records"
done

duplicate_count="$(awk -F '\t' 'NR > 1 {key=$1 FS $2 FS $3; if (seen[key]++) duplicates++} END {print duplicates+0}' "$all_tmp")"
[[ "$duplicate_count" -eq 0 ]] || {
  echo "Error: merged table contains $duplicate_count duplicate Group/Sample/Genome keys" >&2
  exit 1
}

mv "$all_tmp" "${OUTPUT_DIR}/IS_AS_ES_NS_MAG_TPM_long.tsv"
mv "$summary_tmp" "${OUTPUT_DIR}/IS_AS_ES_NS_MAG_TPM_summary.tsv"

printf '[%s] Merge complete: %s\n' "$(date '+%F %T')" \
  "${OUTPUT_DIR}/IS_AS_ES_NS_MAG_TPM_long.tsv"
printf '[%s] Summary: %s\n' "$(date '+%F %T')" \
  "${OUTPUT_DIR}/IS_AS_ES_NS_MAG_TPM_summary.tsv"
