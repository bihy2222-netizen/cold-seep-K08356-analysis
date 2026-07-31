#!/usr/bin/env bash
set -euo pipefail

# Extract MAG proteins by direct hmmsearch against a target HMM model.
#
# Example matching the original dmdA-style command:
#   HMM_FILE=dmdA.hmm \
#   PROTEINS=selected_MAG_proteins.faa \
#   PREFIX=dmdA \
#   CPU=200 \
#   EVALUE=1e-50 \
#   bash 05_extract_MAG_proteins_by_hmmsearch.sh
#
# Ordered workflow:
#   01. Run hmmsearch against the selected MAG protein FASTA.
#   02. Parse --tblout and keep non-comment hits with full-sequence E-value <= EVALUE.
#   03. Extract matched protein sequences from the input FASTA.
#
# Notes:
#   - In HMMER --tblout, column 1 is target sequence ID.
#   - Column 5 is the full-sequence E-value.
#   - Use this as a direct HMM evidence route, complementary to strict KOfam
#     extraction from first-column "*" rows.

HMM_FILE="${HMM_FILE:-dmdA.hmm}"
PROTEINS="${PROTEINS:-selected_MAG_proteins.faa}"
PREFIX="${PREFIX:-dmdA}"
CPU="${CPU:-200}"
EVALUE="${EVALUE:-1e-50}"
HMMSEARCH_BIN="${HMMSEARCH_BIN:-hmmsearch}"
SEQKIT_BIN="${SEQKIT_BIN:-seqkit}"
OUT_DIR="${OUT_DIR:-.}"

mkdir -p "$OUT_DIR"

TBL_OUT="${OUT_DIR}/${PREFIX}_hits.tbl"
HMM_OUT="${OUT_DIR}/${PREFIX}_hmm.out"
IDS_OUT="${OUT_DIR}/${PREFIX}.ids"
FAA_OUT="${OUT_DIR}/${PREFIX}_candidates.faa"

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

require_program "$HMMSEARCH_BIN"
require_program "$SEQKIT_BIN"

[[ -s "$HMM_FILE" ]] || { echo "Error: HMM file not found: $HMM_FILE" >&2; exit 1; }
[[ -s "$PROTEINS" ]] || { echo "Error: protein FASTA not found: $PROTEINS" >&2; exit 1; }

echo "[$(date '+%F %T')] 01. hmmsearch: $PREFIX"
"$HMMSEARCH_BIN" \
  --cpu "$CPU" \
  -E "$EVALUE" \
  --tblout "$TBL_OUT" \
  "$HMM_FILE" \
  "$PROTEINS" \
  > "$HMM_OUT"

echo "[$(date '+%F %T')] 02. Parse tblout by full-sequence E-value <= $EVALUE"
awk -v evalue="$EVALUE" '
  $1 !~ /^#/ && ($5 + 0) <= (evalue + 0) {
    print $1
  }
' "$TBL_OUT" | sort -u > "$IDS_OUT"

echo "[$(date '+%F %T')] 03. Extract protein sequences"
if [[ ! -s "$IDS_OUT" ]]; then
  : > "$FAA_OUT"
  echo "No hits passed the threshold. Empty FASTA written: $FAA_OUT"
else
  "$SEQKIT_BIN" grep -f "$IDS_OUT" "$PROTEINS" > "$FAA_OUT"
fi

printf 'HMM\t%s\n' "$HMM_FILE"
printf 'Protein FASTA\t%s\n' "$PROTEINS"
printf 'E-value cutoff\t%s\n' "$EVALUE"
printf 'Hit IDs\t%s\t%s\n' "$IDS_OUT" "$(wc -l < "$IDS_OUT" | tr -d ' ')"
printf 'Extracted FASTA\t%s\t%s sequences\n' "$FAA_OUT" "$(count_fasta "$FAA_OUT")"
