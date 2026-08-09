#!/usr/bin/env bash
set -euo pipefail

FULL_FAA=${1:?Usage: $0 FULL_808MAG_FAA MODEL_DIR WORKDIR [THREADS]}
MODEL_DIR=${2:?Usage: $0 FULL_808MAG_FAA MODEL_DIR WORKDIR [THREADS]}
WORKDIR=${3:?Usage: $0 FULL_808MAG_FAA MODEL_DIR WORKDIR [THREADS]}
THREADS=${4:-32}
HMMSEARCH=${HMMSEARCH:-hmmsearch}

mkdir -p "$WORKDIR"/{00_input,01_qc,02_rescore58,03_full808_T100,04_T100_union_rescore,05_tables,logs}
seqkit stats -a -T "$FULL_FAA" > "$WORKDIR/01_qc/full808_proteins.seqkit_stats.tsv"
sha256sum "$FULL_FAA" > "$WORKDIR/01_qc/full808_input.sha256"

awk '/^>/{n++; sub(/^>/,""); split($1,a,"|"); mags[a[1]]=1}
  END{for(m in mags)u++; print "protein_count\t" n; print "unique_MAG_count\t" u}' \
  "$FULL_FAA" > "$WORKDIR/01_qc/full808_counts.tsv"

run_models() {
  local input=$1 outdir=$2 prefix=$3 mode=$4
  local name model
  for spec in \
    "IdrA iriA_new.hmm" \
    "AioA aioA.hmm" \
    "combined combined_iriA_aioA.hmm"; do
    read -r name model <<< "$spec"
    if [[ $mode == max ]]; then
      "$HMMSEARCH" --max --cpu "$THREADS" -T 0 --domT 0 \
        --tblout "$outdir/$prefix.$name.tbl" \
        --domtblout "$outdir/$prefix.$name.domtbl" \
        "$MODEL_DIR/$model" "$input" > "$outdir/$prefix.$name.out" \
        2> "$WORKDIR/logs/$prefix.$name.stderr.log"
    else
      "$HMMSEARCH" --cpu "$THREADS" -T 100 --domT 20 \
        --tblout "$outdir/$prefix.$name.T100.tbl" \
        --domtblout "$outdir/$prefix.$name.T100.domtbl" \
        "$MODEL_DIR/$model" "$input" > "$outdir/$prefix.$name.T100.out" \
        2> "$WORKDIR/logs/$prefix.$name.T100.stderr.log"
    fi
  done
}

# Optional urgent stage: place an existing 58-protein union at this path.
if [[ -s "$WORKDIR/00_input/T640_union58.faa" ]]; then
  run_models "$WORKDIR/00_input/T640_union58.faa" "$WORKDIR/02_rescore58" 58_union max
fi

run_models "$FULL_FAA" "$WORKDIR/03_full808_T100" full808 discovery
awk '!/^#/{print $1}' "$WORKDIR"/03_full808_T100/*.tbl | sort -u \
  > "$WORKDIR/00_input/T100_union_prefixed_ids.txt"
seqkit grep -f "$WORKDIR/00_input/T100_union_prefixed_ids.txt" "$FULL_FAA" \
  -o "$WORKDIR/00_input/T100_union.faa"
run_models "$WORKDIR/00_input/T100_union.faa" "$WORKDIR/04_T100_union_rescore" T100_union max
