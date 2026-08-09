#!/usr/bin/env bash
set -euo pipefail

# Standardized evidence workflow:
# 1) author combined HMM, hard threshold -T 640
# 2) separate combined/IdrA/AioA scores for K08356 candidates
# 3) +/-10 ORF extraction with explicit edge censoring
# 4) literature-validated IdrB/AioB-related + P-like reference search
# 5) independent synteny_status and non-circular synteny score

ROOT="${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep}"
ANN="${ANN:-${ROOT}/illu/binning/K08356_four_group_annotation}"
HMM="${HMM:-${ROOT}/Neighborhood_Analyses/hmm}"
REF="${REF:-${ROOT}/Idr_reference_tree_assessment/results}"
THREADS="${THREADS:-20}"
STAMP="${STAMP:-$(date '+%Y%m%d_%H%M%S')}"
OUT="${OUT:-${ANN}/standardized_DMSOR_evidence_${STAMP}}"
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
DIAMOND="${DIAMOND:-/home/ps/anaconda3/envs/diamond/bin/diamond}"

mkdir -p "$OUT"/{00_input,01_hmm,02_neighborhood,03_reference_search,04_tables,logs}

CAND="${ANN}/all_groups.K08356.faa"
ALL="$OUT/00_input/all_groups.all_derep_MAGs.faa"
BREF="${REF}/IdrB_AioB.dedup.faa"
PREF="${REF}/IdrP_like.dedup.faa"
JREF="$OUT/00_input/validated_IdrB_AioB_P_like_joint_reference.faa"

for f in "$CAND" "$HMM/combined_iriA_aioA.hmm" "$HMM/iriA_new.hmm" "$HMM/aioA.hmm" "$BREF" "$PREF"; do
  test -s "$f" || { echo "missing/empty: $f" >&2; exit 1; }
done
for x in hmmsearch python3; do command -v "$x" >/dev/null || { echo "missing: $x" >&2; exit 1; }; done
test -x "$DIAMOND" || { echo "missing DIAMOND: $DIAMOND" >&2; exit 1; }

cat "$ANN"/{IS,AS,ES,NS}/01_prodigal/*.all_derep_MAGs.faa > "$ALL"
cat "$BREF" "$PREF" > "$JREF"

{
  date
  sha256sum "$CAND" "$ALL" "$HMM/combined_iriA_aioA.hmm" "$HMM/iriA_new.hmm" "$HMM/aioA.hmm" "$BREF" "$PREF"
  printf 'K08356_candidates\t'; grep -c '^>' "$CAND"
  printf 'all_derep_MAG_proteins\t'; grep -c '^>' "$ALL"
} > "$OUT/00_input/manifest_checksums.tsv"

hmmsearch --cpu "$THREADS" -T 640 --tblout "$OUT/01_hmm/all_derep_combined_T640.tbl" \
  "$HMM/combined_iriA_aioA.hmm" "$ALL" > "$OUT/01_hmm/all_derep_combined_T640.out"
hmmsearch --cpu "$THREADS" --tblout "$OUT/01_hmm/K08356_combined.tbl" \
  "$HMM/combined_iriA_aioA.hmm" "$CAND" > "$OUT/01_hmm/K08356_combined.out"
hmmsearch --cpu "$THREADS" --tblout "$OUT/01_hmm/K08356_IdrA.tbl" \
  "$HMM/iriA_new.hmm" "$CAND" > "$OUT/01_hmm/K08356_IdrA.out"
hmmsearch --cpu "$THREADS" --tblout "$OUT/01_hmm/K08356_AioA.tbl" \
  "$HMM/aioA.hmm" "$CAND" > "$OUT/01_hmm/K08356_AioA.out"

python3 "$HERE/extract_neighborhoods.py" --annotation-root "$ANN" --targets "$CAND" \
  --out-long "$OUT/02_neighborhood/neighborhood_long.tsv" \
  --out-summary "$OUT/02_neighborhood/neighborhood_observability.tsv"

awk 'BEGIN{FS=OFS="\t"} NR==FNR{if(/^>/){id=substr($1,2);keep=0}else if(keep){};next}' /dev/null /dev/null >/dev/null
python3 - "$OUT/02_neighborhood/neighborhood_long.tsv" "$ALL" "$OUT/02_neighborhood/neighborhood_proteins.faa" <<'PY'
import csv,sys
want={r["neighbor_raw_id"] for r in csv.DictReader(open(sys.argv[1]),delimiter="\t")}
keep=False
with open(sys.argv[2],errors="ignore") as i,open(sys.argv[3],"w") as o:
 for line in i:
  if line.startswith(">"): keep=line[1:].split()[0] in want
  if keep:o.write(line)
PY

"$DIAMOND" makedb --in "$JREF" -d "$OUT/03_reference_search/validated_joint_reference" \
  > "$OUT/logs/diamond_makedb.log" 2>&1
"$DIAMOND" blastp --query "$OUT/02_neighborhood/neighborhood_proteins.faa" \
  --db "$OUT/03_reference_search/validated_joint_reference.dmnd" \
  --out "$OUT/03_reference_search/joint_hits.raw.tsv" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
  --more-sensitive --evalue 1e-5 --max-target-seqs 25 --threads "$THREADS" \
  > "$OUT/logs/diamond_blastp.log" 2>&1

python3 "$HERE/summarize_synteny.py" \
  --neighborhood-long "$OUT/02_neighborhood/neighborhood_long.tsv" \
  --neighborhood-summary "$OUT/02_neighborhood/neighborhood_observability.tsv" \
  --diamond "$OUT/03_reference_search/joint_hits.raw.tsv" \
  --out-hits "$OUT/03_reference_search/joint_hits.filtered.tsv" \
  --out-summary "$OUT/04_tables/synteny_evidence.tsv"

printf 'output\t%s\n' "$OUT"
printf 'T640_hits\t'; awk '!/^#/&&NF{n++}END{print n+0}' "$OUT/01_hmm/all_derep_combined_T640.tbl"
printf 'synteny_rows\t'; awk 'END{print NR-1}' "$OUT/04_tables/synteny_evidence.tsv"
