#!/usr/bin/env bash
set -euo pipefail

# Unified MAFFT -> trimAl -> IQ-TREE workflow for the large MopB tree
# and the smaller AioA/IdrA tree.
#
# Usage:
#   bash run_MopB_AioA_IdrA_phylogeny.sh
#   bash run_MopB_AioA_IdrA_phylogeny.sh MopB_all.fasta AioA_IdrA_all.fasta
#
# Input files default to the same directory as this script.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MAFFT_BIN="${MAFFT_BIN:-mafft}"
TRIMAL_BIN="${TRIMAL_BIN:-trimal}"
IQTREE_BIN="${IQTREE_BIN:-iqtree}"

MOPB_INPUT="${1:-MopB_all.fasta}"
AIOA_IDRA_INPUT="${2:-AioA_IdrA_all.fasta}"

for program in "$MAFFT_BIN" "$TRIMAL_BIN" "$IQTREE_BIN"; do
  if ! command -v "$program" >/dev/null 2>&1; then
    echo "Error: required program not found: $program" >&2
    exit 1
  fi
done

for input in "$MOPB_INPUT" "$AIOA_IDRA_INPUT"; do
  if [[ ! -s "$input" ]]; then
    echo "Error: input FASTA not found or empty: $SCRIPT_DIR/$input" >&2
    exit 1
  fi
done

prepare_alignment() {
  local input="$1"
  local prefix="$2"
  local aligned="${prefix}.aligned.fasta"
  local trimmed="${prefix}.trimmed.fasta"

  echo "[$(date '+%F %T')] MAFFT: $input"
  "$MAFFT_BIN" --auto --quiet --anysymbol "$input" > "$aligned"

  echo "[$(date '+%F %T')] trimAl: $aligned"
  "$TRIMAL_BIN" \
    -in "$aligned" \
    -out "$trimmed" \
    -gt 0.8 \
    -st 0.001 \
    -cons 60
}

launch_iqtree() {
  local prefix="$1"
  local trimmed="${prefix}.trimmed.fasta"
  local log="${prefix}.iqtree.log"
  local pid_file="${prefix}.iqtree.pid"

  echo "[$(date '+%F %T')] Launching IQ-TREE: $trimmed"
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
  echo "Log: $SCRIPT_DIR/$log"
}

# Complete both preprocessing workflows before launching the two tree jobs.
prepare_alignment "$MOPB_INPUT" "MopB_all"
prepare_alignment "$AIOA_IDRA_INPUT" "AioA_IdrA_all"

launch_iqtree "MopB_all"
launch_iqtree "AioA_IdrA_all"

echo "Both IQ-TREE jobs are running in the background."
echo "Check them with: ps -p \$(cat MopB_all.iqtree.pid) -p \$(cat AioA_IdrA_all.iqtree.pid)"
echo "Follow logs with: tail -f MopB_all.iqtree.log AioA_IdrA_all.iqtree.log"
