#!/usr/bin/env bash
set -euo pipefail

INPUT_FAA=${1:?Usage: $0 CORE49_FAA WORKDIR [THREADS]}
WORKDIR=${2:?Usage: $0 CORE49_FAA WORKDIR [THREADS]}
THREADS=${3:-20}

MAFFT_BIN=${MAFFT_BIN:-mafft}
TRIMAL_BIN=${TRIMAL_BIN:-trimal}
MMSEQS_BIN=${MMSEQS_BIN:-mmseqs}
NEEDLEALL_BIN=${NEEDLEALL_BIN:-needleall}

mkdir -p "$WORKDIR"/{00_input,01_qc,02_mafft,03_global_identity,04_local_identity,logs}
ln -sfn "$(realpath "$INPUT_FAA")" "$WORKDIR/00_input/K08356_49_sequences.faa"

seqkit stats -a -T "$INPUT_FAA" > "$WORKDIR/01_qc/seqkit_stats.tsv"
seqkit fx2tab -n -l "$INPUT_FAA" > "$WORKDIR/01_qc/sequence_lengths.tsv"
sha256sum "$INPUT_FAA" > "$WORKDIR/01_qc/input.sha256"

"$MAFFT_BIN" --localpair --maxiterate 1000 --thread "$THREADS" "$INPUT_FAA" \
  > "$WORKDIR/02_mafft/K08356_49_sequences.mafft_linsi.faa" \
  2> "$WORKDIR/logs/mafft_linsi.log"

"$TRIMAL_BIN" \
  -in "$WORKDIR/02_mafft/K08356_49_sequences.mafft_linsi.faa" \
  -out "$WORKDIR/02_mafft/K08356_49_sequences.mafft_linsi.trimmed.faa" \
  -automated1 > "$WORKDIR/logs/trimal.log" 2>&1

"$MMSEQS_BIN" easy-search "$INPUT_FAA" "$INPUT_FAA" \
  "$WORKDIR/04_local_identity/mmseqs_all_vs_all.tsv" \
  "$WORKDIR/04_local_identity/mmseqs_tmp" \
  --threads "$THREADS" --max-seqs 1000 -s 7.5 \
  --format-output 'query,target,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits' \
  > "$WORKDIR/logs/mmseqs.log" 2>&1

"$NEEDLEALL_BIN" -asequence "$INPUT_FAA" -bsequence "$INPUT_FAA" \
  -gapopen 10 -gapextend 0.5 -endweight Y -endopen 10 -endextend 0.5 \
  -aformat3 pair -outfile "$WORKDIR/03_global_identity/needleall_all_vs_all.pair" \
  > "$WORKDIR/logs/needleall.stdout.log" 2> "$WORKDIR/logs/needleall.stderr.log"
