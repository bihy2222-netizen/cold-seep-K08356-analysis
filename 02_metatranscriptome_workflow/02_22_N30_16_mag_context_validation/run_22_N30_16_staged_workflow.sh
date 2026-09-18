#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep}
DONG=${DONG:-$ROOT/public_validation/Dong_NewbornSeep_PRJNA1313825_20260803}
PILOT=${PILOT:-$DONG/pilot4_v2}
SAMPLE_DIR="$PILOT/22_N30_16"
WORK="$SAMPLE_DIR/10_maggie_mag_context_validation"
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$WORK"/{00_audit,01_candidate_reconciliation,02_neighborhood,03_contig_bin_assignment,04_personalized_references,05_dna_presence_control_mapping,06_per_base_depth,07_figures,08_igv,09_final_tables}

case "${1:-audit}" in
  audit)
    python3 "$HERE/scripts/01_audit_22_N30_16_inputs.py" --pilot-root "$PILOT" --outdir "$WORK/00_audit"
    echo "Audit complete. Resolve every MISSING/AMBIGUOUS file class before reconciliation."
    ;;
  reconcile)
    : "${COMBINED_TBL:?}" "${COMBINED_T640_TBL:?}" "${COMBINED_DOMTBL:?}" "${IRIA_TBL:?}" "${AIOA_TBL:?}"
    : "${FAA:?}" "${FNA:?}" "${GFF:?}" "${CONTIGS:?}"
    python3 "$HERE/scripts/02_reconcile_22_N30_16_candidates.py" \
      --combined "$COMBINED_TBL" --combined-strict "$COMBINED_T640_TBL" --combined-domtbl "$COMBINED_DOMTBL" \
      --iria "$IRIA_TBL" --aioa "$AIOA_TBL" --faa "$FAA" --fna "$FNA" --gff "$GFF" --contigs "$CONTIGS" \
      --outdir "$WORK/01_candidate_reconciliation"
    echo "Reconciliation complete. Perform tree review and validated IdrB/P-like annotation next."
    ;;
  annotate-neighborhood)
    : "${IDRB_REF:?}" "${PLIKE_REF:?}"
    NEIGHBORHOOD_FAA=${NEIGHBORHOOD_FAA:-$WORK/01_candidate_reconciliation/22_N30_16_candidate_neighborhood_proteins.faa}
    command -v diamond >/dev/null
    cat "$IDRB_REF" "$PLIKE_REF" > "$WORK/02_neighborhood/validated_IdrB_P_like_refs.faa"
    diamond makedb --in "$WORK/02_neighborhood/validated_IdrB_P_like_refs.faa" -d "$WORK/02_neighborhood/validated_IdrB_P_like_refs"
    diamond blastp --query "$NEIGHBORHOOD_FAA" --db "$WORK/02_neighborhood/validated_IdrB_P_like_refs.dmnd" \
      --out "$WORK/02_neighborhood/validated_neighbor_hits.raw.tsv" \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --more-sensitive --evalue 1e-5 --max-target-seqs 25 --threads "${THREADS:-16}"
    python3 "$HERE/scripts/annotate_neighborhood_with_validated_refs.py" \
      --neighborhood "$WORK/01_candidate_reconciliation/22_N30_16_candidate_neighborhood.tsv" \
      --diamond "$WORK/02_neighborhood/validated_neighbor_hits.raw.tsv" \
      --out-neighborhood "$WORK/02_neighborhood/22_N30_16_candidate_neighborhood.reviewed.tsv" \
      --out-hits "$WORK/02_neighborhood/22_N30_16_validated_neighbor_hits.tsv"
    python3 "$HERE/scripts/plot_candidate_neighborhoods.py" \
      --neighborhood "$WORK/02_neighborhood/22_N30_16_candidate_neighborhood.reviewed.tsv" \
      --outdir "$WORK/07_figures/neighborhoods"
    ;;
  finalize)
    : "${FAA:?}" "${FNA:?}" "${CONTIGS:?}" "${BIN_DIR:?}"
    REVIEW=${REVIEW:-$WORK/01_candidate_reconciliation/22_N30_16_candidate_review.tsv}
    NEIGHBORHOOD=${NEIGHBORHOOD:-$WORK/02_neighborhood/22_N30_16_candidate_neighborhood.reviewed.tsv}
    EXTRA=(); [[ -n "${MAG_METADATA:-}" ]] && EXTRA+=(--mag-metadata "$MAG_METADATA")
    python3 "$HERE/scripts/03_finalize_22_N30_16_candidates.py" --review "$REVIEW" --neighborhood "$NEIGHBORHOOD" \
      --faa "$FAA" --fna "$FNA" --contigs "$CONTIGS" --bin-dir "$BIN_DIR" \
      --outdir "$WORK/04_personalized_references" "${EXTRA[@]}"
    ;;
  dna-map)
    : "${CLEAN_R1:?}" "${CLEAN_R2:?}" "${PRODIGAL_GFF:?}" "${WHOLE_ASSEMBLY:?}"
    CLEAN_R1="$CLEAN_R1" CLEAN_R2="$CLEAN_R2" PRODIGAL_GFF="$PRODIGAL_GFF" WHOLE_ASSEMBLY="$WHOLE_ASSEMBLY" \
      PERSONALIZED_DIR="$WORK/04_personalized_references" WORK="$WORK" \
      bash "$HERE/scripts/04_run_22_N30_16_DNA_positive_control.sh"
    ;;
  *) echo "Usage: $0 {audit|reconcile|annotate-neighborhood|finalize|dna-map}" >&2; exit 2 ;;
esac
