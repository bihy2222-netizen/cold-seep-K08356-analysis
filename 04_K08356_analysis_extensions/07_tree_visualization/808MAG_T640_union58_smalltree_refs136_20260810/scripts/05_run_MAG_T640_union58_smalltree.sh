#!/usr/bin/env bash
set -euo pipefail

OUT="/home/ps/ps1/data/bihongyu/cold_seep/HMM_supplement/808MAG_T640_union58_smalltree_refs136_20260810"
CAND="/home/ps/ps1/data/bihongyu/cold_seep/HMM_supplement/808MAG_IdrA_rescore_20260809similarity/00_input/T640_union58.faa"
REF="/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/K08356_four_group_annotation/strict_kofam_best/tree_small_AioA_IdrA_Unknown_group_derep49_refs136/reference_AioA_IdrA_Unknown_136_from_MAG_idratree4_small.faa"
MAFFT="/home/ps/anaconda3/envs/mafft/bin/mafft"
TRIMAL="/home/ps/anaconda3/envs/trimal/bin/trimal"
IQTREE="/home/ps/anaconda3/envs/iqtree/bin/iqtree"
PREFIX="MAG_T640_union58_plus_refs136"

mkdir -p "$OUT"/{00_input,01_alignment,02_tree,03_audit,logs}
exec > >(tee -a "$OUT/logs/workflow.log") 2>&1

echo "[$(date)] workflow started"
cp -f "$CAND" "$OUT/00_input/T640_union58.faa"
cp -f "$REF" "$OUT/00_input/smalltree_refs136.faa"
cat "$OUT/00_input/T640_union58.faa" "$OUT/00_input/smalltree_refs136.faa" \
  > "$OUT/00_input/${PREFIX}.faa"

grep -c '^>' "$OUT/00_input/T640_union58.faa" > "$OUT/03_audit/candidate_sequence_count.txt"
grep -c '^>' "$OUT/00_input/smalltree_refs136.faa" > "$OUT/03_audit/reference_sequence_count.txt"
grep -c '^>' "$OUT/00_input/${PREFIX}.faa" > "$OUT/03_audit/combined_sequence_count.txt"
sha256sum "$OUT/00_input/"*.faa > "$OUT/03_audit/input_sha256.tsv"

awk '/^>/{id=substr($1,2); n[id]++} END{for(id in n) if(n[id]>1) print id, n[id]}' \
  "$OUT/00_input/${PREFIX}.faa" > "$OUT/03_audit/duplicate_first_token_ids.tsv"
if [[ -s "$OUT/03_audit/duplicate_first_token_ids.tsv" ]]; then
  echo "Duplicate first-token FASTA IDs detected; stopping before alignment." >&2
  exit 2
fi

"$MAFFT" --auto --thread 16 "$OUT/00_input/${PREFIX}.faa" \
  > "$OUT/01_alignment/${PREFIX}.aligned.fasta" \
  2> "$OUT/logs/mafft.stderr.log"

"$TRIMAL" -in "$OUT/01_alignment/${PREFIX}.aligned.fasta" \
  -out "$OUT/01_alignment/${PREFIX}.trimmed.fasta" -automated1 \
  > "$OUT/logs/trimal.stdout.log" 2> "$OUT/logs/trimal.stderr.log"

cd "$OUT/02_tree"
"$IQTREE" -s "$OUT/01_alignment/${PREFIX}.trimmed.fasta" \
  -m MFP -B 1000 --alrt 1000 -T AUTO -pre "$PREFIX"

{
  echo -e "candidate_sequences\t$(cat "$OUT/03_audit/candidate_sequence_count.txt")"
  echo -e "reference_sequences\t$(cat "$OUT/03_audit/reference_sequence_count.txt")"
  echo -e "combined_sequences\t$(cat "$OUT/03_audit/combined_sequence_count.txt")"
  echo -e "treefile\t$OUT/02_tree/${PREFIX}.treefile"
  echo -e "completed_at\t$(date -Iseconds)"
} > "$OUT/03_audit/run_summary.tsv"

echo "[$(date)] workflow completed"
