#!/usr/bin/env bash
set -euo pipefail

# Strictly extract group-dereplicated MAG K08356/AioA proteins from KOfamScan.
#
# This follows the original strict idea:
#   01. Merge each group's *.kofamscan.txt files.
#   02. Keep only reliable KOfamScan rows where the first column is "*".
#   03. Keep K08356 rows from the KO column.
#   04. Split "sample_binN-gene_id" into bin_id and gene_id for the summary table.
#   05. Extract matching protein sequences from the corresponding *.cds.faa files.
#
# Server default:
#   /home/ps/ps1/data/bihongyu/cold_seep/illu/binning

ROOT_DIR="${ROOT_DIR:-/home/ps/ps1/data/bihongyu/cold_seep/illu/binning}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/K08356_four_group_annotation/strict_kofam_best}"
TARGET_KO="${TARGET_KO:-K08356}"
SEQKIT_BIN="${SEQKIT_BIN:-seqkit}"

GROUPS=(IS AS ES NS)

group_dir() {
  local group="$1"
  case "$group" in
    IS) printf '%s/IS/IS_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    AS) printf '%s/AS/AS_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    ES) printf '%s/ES/ES_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    NS) printf '%s/NS/NS_50_10_dRep/dereplicated_genomes\n' "$ROOT_DIR" ;;
    *) echo "Unknown group: $group" >&2; return 1 ;;
  esac
}

require_program() {
  local program="$1"
  if ! command -v "$program" >/dev/null 2>&1; then
    echo "Error: required program not found: $program" >&2
    exit 1
  fi
}

count_fasta() {
  local fasta="$1"
  if [[ -s "$fasta" ]]; then
    grep -c '^>' "$fasta"
  else
    printf '0\n'
  fi
}

require_program "$SEQKIT_BIN"
mkdir -p "$OUT_DIR"

all_hits="${OUT_DIR}/all_groups.${TARGET_KO}_hits.tsv"
all_faa="${OUT_DIR}/all_groups.${TARGET_KO}.faa"
: > "$all_hits"
: > "$all_faa"

for group in "${GROUPS[@]}"; do
  d="$(group_dir "$group")"
  group_out="${OUT_DIR}/${group}"
  mkdir -p "$group_out"

  [[ -d "$d" ]] || { echo "Error: group directory not found: $d" >&2; exit 1; }
  compgen -G "${d}/*.kofamscan.txt" >/dev/null || {
    echo "Error: no *.kofamscan.txt files found in $d" >&2
    exit 1
  }
  compgen -G "${d}/*.cds.faa" >/dev/null || {
    echo "Error: no *.cds.faa files found in $d" >&2
    exit 1
  }

  cat "${d}"/*.kofamscan.txt > "${group_out}/all_bin.kofamscan.txt"
  awk '$1 == "*"' \
    "${group_out}/all_bin.kofamscan.txt" \
    > "${group_out}/all_bin.kofamscan_filtered.txt"

  awk -F'\t' 'BEGIN{OFS="\t"}
    {
      full_id = $2
      bin_part = full_id
      gene_part = ""
      if (match(full_id, /_bin[0-9]+-/)) {
        bin_part = substr(full_id, 1, RSTART + RLENGTH - 2)
        gene_part = substr(full_id, RSTART + RLENGTH)
      }
      print $0, bin_part, gene_part
    }' \
    "${group_out}/all_bin.kofamscan_filtered.txt" \
    > "${group_out}/all_bin.kofamscan_filtered_split.txt"

  awk -F'\t' -v group="$group" -v ko="$TARGET_KO" 'BEGIN{
      OFS="\t"
      print "group","full_gene_id","bin_id","gene_id","KO","threshold","score","evalue","definition"
    }
    $1 == "*" && $3 == ko {
      full_id = $2
      bin_part = full_id
      gene_part = ""
      if (match(full_id, /_bin[0-9]+-/)) {
        bin_part = substr(full_id, 1, RSTART + RLENGTH - 2)
        gene_part = substr(full_id, RSTART + RLENGTH)
      }
      print group, full_id, bin_part, gene_part, $3, $4, $5, $6, $7
    }' \
    "${group_out}/all_bin.kofamscan.txt" \
    > "${group_out}/${group}.${TARGET_KO}_hits.tsv"

  awk -F'\t' 'NR > 1 {print $2}' \
    "${group_out}/${group}.${TARGET_KO}_hits.tsv" \
    | sort -u \
    > "${group_out}/${group}.${TARGET_KO}.protein_ids.txt"

  "$SEQKIT_BIN" grep \
    -f "${group_out}/${group}.${TARGET_KO}.protein_ids.txt" \
    "${d}"/*.cds.faa \
    | awk -v group="$group" '/^>/{sub(/^>/, ">" group "|")} {print}' \
    > "${group_out}/${group}.${TARGET_KO}.faa"

  if [[ ! -s "$all_hits" ]]; then
    cat "${group_out}/${group}.${TARGET_KO}_hits.tsv" > "$all_hits"
  else
    tail -n +2 "${group_out}/${group}.${TARGET_KO}_hits.tsv" >> "$all_hits"
  fi
  cat "${group_out}/${group}.${TARGET_KO}.faa" >> "$all_faa"

  printf '%s\thits=%s\tfaa=%s\n' \
    "$group" \
    "$(tail -n +2 "${group_out}/${group}.${TARGET_KO}_hits.tsv" | wc -l | tr -d ' ')" \
    "$(count_fasta "${group_out}/${group}.${TARGET_KO}.faa")"
done

printf 'ALL\thits=%s\tfaa=%s\n' \
  "$(tail -n +2 "$all_hits" | wc -l | tr -d ' ')" \
  "$(count_fasta "$all_faa")"

echo "Result directory: $OUT_DIR"
echo "Combined hit table: $all_hits"
echo "Combined protein FASTA: $all_faa"
