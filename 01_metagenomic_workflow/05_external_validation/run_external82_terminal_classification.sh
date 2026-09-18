#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep}
DONG=${DONG:-$ROOT/public_validation/Dong_NewbornSeep_PRJNA1313825_20260803}
PILOT=${PILOT:-$DONG/pilot4_v2}
OUT=${OUT:-$PILOT/10_external82_terminal_classification_20260814_v3}
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

TREE=${TREE:-$PILOT/08_tree/pilot4_v2_candidates_plus_refs_iqtree.treefile}
TREE_FASTA=${TREE_FASTA:-$PILOT/08_tree/pilot4_v2_candidates_plus_refs.faa}
ANCHORS=${ANCHORS:-$ROOT/HMM_supplement/808MAG_T640_union58_smalltree_refs136_20260810/04_tables/candidate58_four_clade_integrated_evidence.tsv}
IDRB_REF=${IDRB_REF:-$ROOT/Idr_reference_tree_assessment/results/IdrB_AioB.dedup.faa}
PLIKE_REF=${PLIKE_REF:-$ROOT/Idr_reference_tree_assessment/results/IdrP_like.dedup.faa}
CONTIG_BIN=${CONTIG_BIN:-$PILOT/07_reports/pilot4_v2_strict_K08356_contig_to_bin.tsv}
MAG_SUMMARY=${MAG_SUMMARY:-$PILOT/07_reports/pilot4_v2_K08356_positive_bin_summary.tsv}
MAG_TAXONOMY=${MAG_TAXONOMY:-$PILOT/07_reports/pilot4_v2_K08356_positive_MQHQ_MAG_taxonomy.tsv}

mkdir -p "$OUT"
python3 "$HERE/scripts/classify_external82_candidates.py" \
  --pilot-root "$PILOT" \
  --tree "$TREE" \
  --tree-fasta "$TREE_FASTA" \
  --anchor-evidence "$ANCHORS" \
  --idrb-ref "$IDRB_REF" \
  --plike-ref "$PLIKE_REF" \
  --contig-bin "$CONTIG_BIN" \
  --mag-summary "$MAG_SUMMARY" \
  --mag-taxonomy "$MAG_TAXONOMY" \
  --diamond "${DIAMOND:-diamond}" \
  --threads "${THREADS:-16}" \
  --outdir "$OUT"
