#!/usr/bin/env bash
set -euo pipefail

# Standard K08356/AioA protein-tree workflow for group-dereplicated MAGs.
#
# Main purpose:
#   Use strict KOfamScan reliable hits ("*" rows) as the primary K08356/AioA
#   candidate set from IS/AS/ES/NS group-dereplicated MAGs, then build:
#     1. a DMSOR broad tree with broad DMSOR references;
#     2. a small AioA/IdrA/Unknown tree with focused references.
#
# Server defaults match the cold_seep analysis workspace.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

ROOT_DIR="${ROOT_DIR:-/home/ps/ps1/data/bihongyu/cold_seep/illu/binning}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/K08356_four_group_annotation/strict_kofam_best}"
TARGET_KO="${TARGET_KO:-K08356}"

OLD48_FASTA="${OLD48_FASTA:-/home/ps/ps1/data/bihongyu/cold_seep/HMM_supplement/input/48_candidates.faa}"
BROAD_REF_ALL="${BROAD_REF_ALL:-/home/ps/ps1/data/bihongyu/cold_seep/ASproteintree/20260323/all_proteins.faa}"
SMALL_REF_ALL="${SMALL_REF_ALL:-/home/ps/ps1/data/bihongyu/cold_seep/ASproteintree/MAG+idratree4-small-smalltree/all_proteins.faa}"

EXTRACT_SCRIPT="${EXTRACT_SCRIPT:-${SCRIPT_DIR}/04_extract_group_derep_K08356_from_kofamscan.sh}"

SEQKIT_BIN="${SEQKIT_BIN:-seqkit}"
MAFFT_BIN="${MAFFT_BIN:-/home/ps/anaconda3/envs/mafft/bin/mafft}"
TRIMAL_BIN="${TRIMAL_BIN:-/home/ps/anaconda3/envs/trimal/bin/trimal}"
IQTREE_BIN="${IQTREE_BIN:-/home/ps/anaconda3/envs/iqtree/bin/iqtree}"

RUN_EXTRACT="${RUN_EXTRACT:-1}"
RUN_BIG_TREE="${RUN_BIG_TREE:-1}"
RUN_SMALL_TREE="${RUN_SMALL_TREE:-1}"
RUN_IQTREE="${RUN_IQTREE:-1}"

require_file() {
  local file="$1"
  [[ -s "$file" ]] || { echo "Error: required file not found or empty: $file" >&2; exit 1; }
}

require_program() {
  local program="$1"
  [[ -x "$program" ]] || command -v "$program" >/dev/null 2>&1 || {
    echo "Error: required program not found: $program" >&2
    exit 1
  }
}

count_fasta() {
  local fasta="$1"
  if [[ -s "$fasta" ]]; then
    grep -c '^>' "$fasta"
  else
    printf '0\n'
  fi
}

clean_aa() {
  awk '/^>/{print; next} {gsub(/[^ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]/, "X"); print}'
}

prepare_old48_ids() {
  local ids_out="$1"
  grep '^>' "$OLD48_FASTA" | sed 's/^>//; s/ .*//' | sort -u > "$ids_out"
}

run_alignment_and_tree() {
  local workdir="$1"
  local prefix="$2"
  local input="${workdir}/${prefix}.fasta"
  local aligned="${workdir}/${prefix}.aligned.fasta"
  local trimmed="${workdir}/${prefix}.trimmed.fasta"
  local log="${workdir}/${prefix}.iqtree.log"
  local pid_file="${workdir}/${prefix}.iqtree.pid"

  require_file "$input"
  require_program "$MAFFT_BIN"
  require_program "$TRIMAL_BIN"

  echo "[$(date '+%F %T')] MAFFT: $input"
  "$MAFFT_BIN" --auto --quiet --anysymbol "$input" > "$aligned"

  echo "[$(date '+%F %T')] trimAl: $aligned"
  "$TRIMAL_BIN" \
    -in "$aligned" \
    -out "$trimmed" \
    -gt 0.8 \
    -st 0.001 \
    -cons 60

  "$SEQKIT_BIN" stats "$input" "$aligned" "$trimmed" > "${workdir}/${prefix}.seqkit_stats.tsv"
  cat "${workdir}/${prefix}.seqkit_stats.tsv"

  if [[ "$RUN_IQTREE" != "1" ]]; then
    echo "RUN_IQTREE=$RUN_IQTREE; stop after trimAl for $prefix."
    return
  fi

  require_program "$IQTREE_BIN"
  echo "[$(date '+%F %T')] IQ-TREE: $trimmed"
  nohup "$IQTREE_BIN" \
    -s "$trimmed" \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T AUTO \
    > "$log" 2>&1 &

  local pid=$!
  printf '%s\n' "$pid" > "$pid_file"
  echo "$prefix IQ-TREE PID: $pid"
  echo "Log: $log"
}

make_big_tree_input() {
  local mag_faa="$1"
  local old48_ids="$2"
  local workdir="${OUT_DIR}/tree_big_DMSOR_group_derep49_refs261"
  local prefix="K08356_group_derep49_DMSOR_refs261"
  mkdir -p "$workdir"

  "$SEQKIT_BIN" grep -v -f "$old48_ids" "$BROAD_REF_ALL" \
    > "${workdir}/reference_DMSOR_261_from_20260323_all_proteins.faa"

  cat "$mag_faa" "${workdir}/reference_DMSOR_261_from_20260323_all_proteins.faa" \
    | clean_aa \
    > "${workdir}/${prefix}.fasta"

  "$SEQKIT_BIN" stats \
    "$mag_faa" \
    "${workdir}/reference_DMSOR_261_from_20260323_all_proteins.faa" \
    "${workdir}/${prefix}.fasta" \
    > "${workdir}/input.seqkit_stats.tsv"
  cat "${workdir}/input.seqkit_stats.tsv"

  run_alignment_and_tree "$workdir" "$prefix"
}

make_small_tree_input() {
  local mag_faa="$1"
  local old48_ids="$2"
  local workdir="${OUT_DIR}/tree_small_AioA_IdrA_Unknown_group_derep49_refs136"
  local prefix="K08356_group_derep49_AioA_IdrA_Unknown_refs136"
  mkdir -p "$workdir"

  "$SEQKIT_BIN" grep -v -f "$old48_ids" "$SMALL_REF_ALL" \
    > "${workdir}/reference_AioA_IdrA_Unknown_136_from_MAG_idratree4_small.faa"

  cat "$mag_faa" "${workdir}/reference_AioA_IdrA_Unknown_136_from_MAG_idratree4_small.faa" \
    | clean_aa \
    > "${workdir}/${prefix}.fasta"

  "$SEQKIT_BIN" stats \
    "$mag_faa" \
    "${workdir}/reference_AioA_IdrA_Unknown_136_from_MAG_idratree4_small.faa" \
    "${workdir}/${prefix}.fasta" \
    > "${workdir}/input.seqkit_stats.tsv"
  cat "${workdir}/input.seqkit_stats.tsv"

  run_alignment_and_tree "$workdir" "$prefix"
}

main() {
  require_program "$SEQKIT_BIN"
  require_file "$OLD48_FASTA"
  require_file "$BROAD_REF_ALL"
  require_file "$SMALL_REF_ALL"

  mkdir -p "$OUT_DIR"

  if [[ "$RUN_EXTRACT" == "1" ]]; then
    echo "[$(date '+%F %T')] Extract strict group-derep ${TARGET_KO} proteins"
    ROOT_DIR="$ROOT_DIR" OUT_DIR="$OUT_DIR" TARGET_KO="$TARGET_KO" bash "$EXTRACT_SCRIPT"
  fi

  local mag_faa="${OUT_DIR}/all_groups.${TARGET_KO}.faa"
  require_file "$mag_faa"
  echo "Strict group-derep ${TARGET_KO} proteins: $(count_fasta "$mag_faa")"

  local old48_ids="${OUT_DIR}/old48_MAG_ids.txt"
  prepare_old48_ids "$old48_ids"

  if [[ "$RUN_BIG_TREE" == "1" ]]; then
    echo "[$(date '+%F %T')] Build broad DMSOR tree input and launch tree"
    make_big_tree_input "$mag_faa" "$old48_ids"
  fi

  if [[ "$RUN_SMALL_TREE" == "1" ]]; then
    echo "[$(date '+%F %T')] Build small AioA/IdrA/Unknown tree input and launch tree"
    make_small_tree_input "$mag_faa" "$old48_ids"
  fi
}

main "$@"
