#!/usr/bin/env bash
set -euo pipefail

# Run CoverM TPM profiling for IS/AS/ES/NS group-dereplicated MAG sets.
#
# This reproduces the previous overall derep CoverM settings, but changes the
# genome directory and output directory for each group:
#   coverm genome
#     --genome-fasta-extension fasta
#     -m tpm
#     --min-covered-fraction 0.1
#     --min-read-percent-identity 95
#     --min-read-aligned-percent 75
#     --trim-min 0.1
#     --trim-max 0.9
#     -t 100
#
# Server default:
#   /home/ps/ps1/data/bihongyu/cold_seep/illu/binning

ROOT_DIR="${ROOT_DIR:-/home/ps/ps1/data/bihongyu/cold_seep/illu/binning}"
OLD_DEREP_DIR="${OLD_DEREP_DIR:-${ROOT_DIR}/all_50_10_dRep/dereplicated_genomes}"
COVERM_BIN="${COVERM_BIN:-/home/ps/anaconda3/envs/coverm/bin/coverm}"
export PATH="/home/ps/anaconda3/envs/coverm/bin:${PATH}"
THREADS="${THREADS:-100}"

HABITAT_GROUPS=(IS AS ES NS)
SAMPLE_TABLE="${ROOT_DIR}/group_coverm_samples.tsv"

make_sample_table() {
  awk '
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
        print sample "\t" r1 "\t" r2
      }
    }
  ' "${OLD_DEREP_DIR}/coverm1.sh" "${OLD_DEREP_DIR}/coverm2.sh" | sort -u > "$SAMPLE_TABLE"
}

group_genome_dir() {
  local group="$1"
  case "$group" in
    IS) printf '%s/IS/IS_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    AS) printf '%s/AS/AS_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    ES) printf '%s/ES/ES_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    NS) printf '%s/NS/NS_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    *) echo "Unknown group: $group" >&2; return 1 ;;
  esac
}

if [[ ! -x "$COVERM_BIN" ]]; then
  echo "Error: coverm not executable: $COVERM_BIN" >&2
  exit 1
fi

if [[ ! -s "$SAMPLE_TABLE" ]]; then
  make_sample_table
fi

printf '[%s] Sample table: %s (%s samples)\n' \
  "$(date '+%F %T')" \
  "$SAMPLE_TABLE" \
  "$(wc -l < "$SAMPLE_TABLE" | tr -d ' ')"

for group in "${HABITAT_GROUPS[@]}"; do
  genome_dir="$(group_genome_dir "$group")"
  out_dir="${genome_dir}/bin_tpm"
  mkdir -p "$out_dir"

  [[ -d "$genome_dir" ]] || { echo "Error: genome directory not found: $genome_dir" >&2; exit 1; }

  genome_count="$(find "$genome_dir" -maxdepth 1 -type f -name '*.fasta' | wc -l | tr -d ' ')"
  printf '[%s] Group %s: %s genomes -> %s\n' \
    "$(date '+%F %T')" "$group" "$genome_count" "$out_dir"

  while IFS=$'\t' read -r sample r1 r2; do
    [[ -n "$sample" ]] || continue
    out_file="${out_dir}/${sample}_tpm.txt"

    if [[ -s "$out_file" ]]; then
      printf '[%s] SKIP %s %s existing %s\n' \
        "$(date '+%F %T')" "$group" "$sample" "$out_file"
      continue
    fi

    printf '[%s] RUN %s %s\n' "$(date '+%F %T')" "$group" "$sample"
    "$COVERM_BIN" genome \
      --genome-fasta-directory "$genome_dir" \
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
      -o "$out_file"
  done < "$SAMPLE_TABLE"
done

printf '[%s] All group CoverM TPM jobs completed.\n' "$(date '+%F %T')"
