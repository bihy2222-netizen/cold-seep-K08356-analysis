#!/usr/bin/env bash
set -euo pipefail
CAND=${1:?48_candidates.faa}
ALLMAG=${2:?all_MAG_proteins.faa}
HMMDIR=${3:?hmm_dir}
OUT=${4:-result}
CPU=${5:-20}
mkdir -p "$OUT/01_hmm_result" "$OUT/03_tables" "$OUT/04_new_candidates"
for x in hmmsearch python3; do command -v "$x" >/dev/null || { echo "missing $x" >&2; exit 1; }; done
for f in "$CAND" "$ALLMAG" "$HMMDIR/combined_iriA_aioA.hmm" "$HMMDIR/iriA_new.hmm" "$HMMDIR/aioA.hmm"; do [ -s "$f" ] || { echo "missing/empty: $f" >&2; exit 1; }; done
printf '[%s] input counts\n' "$(date)"
grep -c '^>' "$CAND" | awk '{print "48_candidates",$1}'
grep -c '^>' "$ALLMAG" | awk '{print "all_MAG_proteins",$1}'
awk '/^>/{print $1}' "$CAND" | sort | uniq -d > "$OUT/03_tables/duplicate_ids_48.txt"
awk '/^>/{print $1}' "$ALLMAG" | sort | uniq -d > "$OUT/03_tables/duplicate_ids_all_MAG.txt"
printf '[%s] hmmsearch combined 48\n' "$(date)"
hmmsearch --cpu "$CPU" --tblout "$OUT/01_hmm_result/combined_48.tbl" "$HMMDIR/combined_iriA_aioA.hmm" "$CAND" > "$OUT/01_hmm_result/combined_48.out"
printf '[%s] hmmsearch idrA 48\n' "$(date)"
hmmsearch --cpu "$CPU" --tblout "$OUT/01_hmm_result/idrA_48.tbl" "$HMMDIR/iriA_new.hmm" "$CAND" > "$OUT/01_hmm_result/idrA_48.out"
printf '[%s] hmmsearch aioA 48\n' "$(date)"
hmmsearch --cpu "$CPU" --tblout "$OUT/01_hmm_result/aioA_48.tbl" "$HMMDIR/aioA.hmm" "$CAND" > "$OUT/01_hmm_result/aioA_48.out"
printf '[%s] hmmsearch all MAG combined -T 640\n' "$(date)"
hmmsearch --cpu "$CPU" -T 640 --tblout "$OUT/01_hmm_result/all_MAG_combined.tbl" "$HMMDIR/combined_iriA_aioA.hmm" "$ALLMAG" > "$OUT/01_hmm_result/all_MAG_combined.out"
printf '[%s] merge tables\n' "$(date)"
python3 02_merge_hmm_scores.py "$HMMDIR" "$OUT" "$CAND" "$ALLMAG"
printf '[%s] done\n' "$(date)"
