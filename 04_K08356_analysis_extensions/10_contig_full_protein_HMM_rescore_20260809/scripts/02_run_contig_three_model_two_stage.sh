#!/usr/bin/env bash
set -euo pipefail

WORKDIR=${WORKDIR:-/home/ps/ps1/data/bihongyu/cold_seep/HMM_supplement/57sample_contig_IdrA_rescore_20260809}
INPUT_DIR=${INPUT_DIR:-/home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/AA}
MODEL_DIR=${MODEL_DIR:-/home/ps/ps1/data/bihongyu/cold_seep/HMM_supplement/808MAG_IdrA_screen_20260807/01_models}
HMMSEARCH=${HMMSEARCH:-/home/ps/anaconda3/envs/hmmer/bin/hmmsearch}
SEQKIT=${SEQKIT:-/usr/bin/seqkit}
PYTHON=${PYTHON:-/usr/bin/python3}
CPU_DISCOVERY=${CPU_DISCOVERY:-12}
CPU_EXACT=${CPU_EXACT:-8}

mkdir -p "$WORKDIR"/{00_input,01_qc,02_models,03_T100_discovery,04_candidate_union,05_union_exact_rescore,06_tables,07_logs,scripts}
exec > >(tee -a "$WORKDIR/07_logs/master_workflow.log") 2>&1
echo "[$(date -Is)] workflow_start host=$(hostname) pid=$$"

cp "$MODEL_DIR/iriA_new.hmm" "$WORKDIR/02_models/"
cp "$MODEL_DIR/aioA.hmm" "$WORKDIR/02_models/"
cp "$MODEL_DIR/combined_iriA_aioA.hmm" "$WORKDIR/02_models/"
sha256sum "$WORKDIR"/02_models/*.hmm > "$WORKDIR/01_qc/model_sha256.tsv"

"$PYTHON" "$WORKDIR/scripts/01_prepare_prefixed_contig_faa.py" \
  --input-dir "$INPUT_DIR" \
  --output-faa "$WORKDIR/00_input/all_valid_samples.prefixed.faa" \
  --manifest "$WORKDIR/01_qc/sample_FAA_audit.tsv" \
  --excluded-file SY457BB8-12_AA.faa

for spec in "IdrA:iriA_new.hmm" "AioA:aioA.hmm" "combined:combined_iriA_aioA.hmm"; do
  model=${spec%%:*}
  hmm=${spec#*:}
  echo "[$(date -Is)] discovery_start model=$model"
  nice -n 10 "$HMMSEARCH" --cpu "$CPU_DISCOVERY" -T 100 --domT 20 \
    --tblout "$WORKDIR/03_T100_discovery/contig.${model}.T100.tbl" \
    --domtblout "$WORKDIR/03_T100_discovery/contig.${model}.T100.domtbl" \
    "$WORKDIR/02_models/$hmm" "$WORKDIR/00_input/all_valid_samples.prefixed.faa" \
    > "$WORKDIR/03_T100_discovery/contig.${model}.T100.out"
  echo "[$(date -Is)] discovery_done model=$model"
done

awk '!/^#/ && NF {print $1}' "$WORKDIR"/03_T100_discovery/*.T100.tbl \
  | LC_ALL=C sort -u > "$WORKDIR/04_candidate_union/T100_union_ids.txt"
"$SEQKIT" grep -f "$WORKDIR/04_candidate_union/T100_union_ids.txt" \
  "$WORKDIR/00_input/all_valid_samples.prefixed.faa" \
  > "$WORKDIR/04_candidate_union/T100_union.faa"
"$SEQKIT" stats -a -T "$WORKDIR/04_candidate_union/T100_union.faa" \
  > "$WORKDIR/04_candidate_union/T100_union.seqkit_stats.tsv"

for spec in "IdrA:iriA_new.hmm" "AioA:aioA.hmm" "combined:combined_iriA_aioA.hmm"; do
  model=${spec%%:*}
  hmm=${spec#*:}
  echo "[$(date -Is)] exact_rescore_start model=$model"
  nice -n 10 "$HMMSEARCH" --max --cpu "$CPU_EXACT" -T 0 --domT 0 \
    --tblout "$WORKDIR/05_union_exact_rescore/contig_union.${model}.exact.tbl" \
    --domtblout "$WORKDIR/05_union_exact_rescore/contig_union.${model}.exact.domtbl" \
    "$WORKDIR/02_models/$hmm" "$WORKDIR/04_candidate_union/T100_union.faa" \
    > "$WORKDIR/05_union_exact_rescore/contig_union.${model}.exact.out"
  echo "[$(date -Is)] exact_rescore_done model=$model"
done

"$PYTHON" "$WORKDIR/scripts/03_summarize_contig_union_rescore.py" --workdir "$WORKDIR"
echo "[$(date -Is)] workflow_complete"
