#!/usr/bin/env bash
set -euo pipefail

# Standardized local post-processing for the small derep49 AioA/IdrA tree.
# The remote HMM -T 640 run is documented below and its outputs are archived in
# output/remote_T640_20260801. This script regenerates the iTOL-order figures
# and summary tables from the archived HMM tables plus local tree/SVG inputs.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$ROOT/.." && pwd)"
SMALL_TREE_DIR="$PROJECT_ROOT/MAG 补充 idra-tree/分组 derep 的 MAG 蛋白树"

python3 "$ROOT/scripts/make_itol_order_gene_island_outputs.py" \
  --input-svg "$ROOT/input/c0J-IDrhpB2odQ8jRmUo-Q.svg" \
  --treefile "$ROOT/input/AioA_IdrA_Unknown_group_derep49_plus_refs136.trimmed.fasta.treefile" \
  --metadata "$SMALL_TREE_DIR/beautified_small_tree/tip_metadata.tsv" \
  --arrow-manifest "$SMALL_TREE_DIR/small_tree_arrow_plot/small_tree_gene_arrow_manifest.tsv" \
  --hmm-scores "$ROOT/output/remote_T640_20260801/03_tables/HMM_scores_48.tsv" \
  --t640 "$ROOT/output/remote_T640_20260801/03_tables/all_MAG_T640_hits.tsv" \
  --strict-qc "$PROJECT_ROOT/MAG 补充 idra-tree/joint_reference_search_48/strict_synteny_qc_48_reclassified.tsv" \
  --out-dir "$ROOT/output/regenerated_iTOL_order"

cp "$ROOT"/output/regenerated_iTOL_order/*.svg "$ROOT/output/figures/"
cp "$ROOT"/output/regenerated_iTOL_order/*.tsv "$ROOT/output/tables/"

printf 'Regenerated outputs under %s/output\n' "$ROOT"
