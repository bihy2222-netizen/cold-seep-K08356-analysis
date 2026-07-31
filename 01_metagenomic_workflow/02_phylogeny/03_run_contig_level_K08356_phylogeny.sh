#!/usr/bin/env bash
set -euo pipefail

# Reproduce the contig-level K08356/AioA extraction and DMSOR reference tree.
#
# Original server path:
#   /home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/CDS/aio-A_contig_level/tree
#
# Original inputs:
#   all_contig_level_K08356.faa  # extracted contig-level K08356/AioA proteins
#   all_ref_DMSO.faa             # DMSOR reference large-subunit proteins
#
# Ordered workflow:
#   01. Extract K08356 proteins from *.cdhit.kofamscan.best.txt + *.cdhit.cds.faa
#   02. Concatenate per-sample *.K08356.faa into all_contig_level_K08356.faa
#   03. Concatenate all_contig_level_K08356.faa + all_ref_DMSO.faa into all_tree.faa
#   04. Run MAFFT with --auto --quiet --anysymbol
#   05. Run trimAl with -gt 0.8 -st 0.001 -cons 60
#   06. Launch IQ-TREE with -m MFP -B 1000 --alrt 1000 -T AUTO

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

CDS_DIR="${CDS_DIR:-/home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/CDS}"
OUT_DIR="${OUT_DIR:-${CDS_DIR}/aio-A_contig_level}"
TREE_DIR="${TREE_DIR:-${OUT_DIR}/tree}"
TARGET_KO="${TARGET_KO:-K08356}"
REF_FASTA="${REF_FASTA:-${TREE_DIR}/all_ref_DMSO.faa}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
MAFFT_BIN="${MAFFT_BIN:-mafft}"
TRIMAL_BIN="${TRIMAL_BIN:-trimal}"
IQTREE_BIN="${IQTREE_BIN:-iqtree}"
RUN_IQTREE="${RUN_IQTREE:-1}"

EXTRACT_SCRIPT="${SCRIPT_DIR}/02_extract_contig_level_K08356_faa.py"
CONTIG_FASTA="${TREE_DIR}/all_contig_level_${TARGET_KO}.faa"
TREE_INPUT="${TREE_DIR}/all_tree.faa"
ALIGNED_FASTA="${TREE_DIR}/all_tree.aln.fasta"
TRIMMED_FASTA="${TREE_DIR}/all_tree.trim.fasta"
IQTREE_LOG="${TREE_DIR}/all_tree.iqtree.log"
IQTREE_PID="${TREE_DIR}/all_tree.iqtree.pid"

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

echo "[$(date '+%F %T')] 00. Check inputs and tools"
[[ -d "$CDS_DIR" ]] || { echo "Error: CDS_DIR not found: $CDS_DIR" >&2; exit 1; }
[[ -s "$REF_FASTA" ]] || { echo "Error: REF_FASTA not found: $REF_FASTA" >&2; exit 1; }
[[ -s "$EXTRACT_SCRIPT" ]] || { echo "Error: extraction script not found: $EXTRACT_SCRIPT" >&2; exit 1; }

require_program "$PYTHON_BIN"
require_program "$MAFFT_BIN"
require_program "$TRIMAL_BIN"
if [[ "$RUN_IQTREE" == "1" ]]; then
  require_program "$IQTREE_BIN"
fi

mkdir -p "$OUT_DIR" "$TREE_DIR"

echo "[$(date '+%F %T')] 01. Extract ${TARGET_KO} proteins from KOfamScan best hits"
"$PYTHON_BIN" "$EXTRACT_SCRIPT" \
  --ko "$TARGET_KO" \
  --input-dir "$CDS_DIR" \
  --output-dir "$OUT_DIR"

echo "[$(date '+%F %T')] 02. Merge per-sample ${TARGET_KO} FASTA files"
mapfile -d '' ko_fastas < <(find "$OUT_DIR" -maxdepth 1 -type f -name "*.${TARGET_KO}.faa" -print0 | sort -z)
if [[ "${#ko_fastas[@]}" -eq 0 ]]; then
  echo "Error: no per-sample ${TARGET_KO} FASTA files found in $OUT_DIR" >&2
  exit 1
fi
cat "${ko_fastas[@]}" > "$CONTIG_FASTA"

echo "Contig-level ${TARGET_KO} sequences: $(count_fasta "$CONTIG_FASTA")"

echo "[$(date '+%F %T')] 03. Merge contig-level proteins with DMSOR references"
cat "$CONTIG_FASTA" "$REF_FASTA" \
  | awk '/^>/{print; next} {gsub(/[^ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]/, "X"); print}' \
  > "$TREE_INPUT"

echo "Tree input sequences: $(count_fasta "$TREE_INPUT")"

echo "[$(date '+%F %T')] 04. MAFFT alignment"
"$MAFFT_BIN" --auto --quiet --anysymbol "$TREE_INPUT" > "$ALIGNED_FASTA"

echo "[$(date '+%F %T')] 05. trimAl filtering"
"$TRIMAL_BIN" \
  -in "$ALIGNED_FASTA" \
  -out "$TRIMMED_FASTA" \
  -gt 0.8 \
  -st 0.001 \
  -cons 60

echo "Trimmed alignment sequences: $(count_fasta "$TRIMMED_FASTA")"

if [[ "$RUN_IQTREE" != "1" ]]; then
  echo "RUN_IQTREE=$RUN_IQTREE; stop after trimAl."
  exit 0
fi

echo "[$(date '+%F %T')] 06. Launch IQ-TREE"
cd "$TREE_DIR"
nohup "$IQTREE_BIN" \
  -s "$TRIMMED_FASTA" \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T AUTO \
  > "$IQTREE_LOG" 2>&1 &

pid=$!
printf '%s\n' "$pid" > "$IQTREE_PID"
echo "IQ-TREE PID: $pid"
echo "Log: $IQTREE_LOG"
