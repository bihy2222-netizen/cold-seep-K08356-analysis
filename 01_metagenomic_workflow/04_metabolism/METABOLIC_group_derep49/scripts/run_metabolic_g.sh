#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_metabolic_g.sh --input-dir DIR --output-dir DIR --log-dir DIR [options]

Required:
  --input-dir DIR
  --output-dir DIR
  --log-dir DIR

Options:
  --threads INT              Default: 100
  --module-cutoff FLOAT      Default: 0.75
  --kofam-db NAME            Default: full
  --conda-env NAME           Default: METABOLIC_v4.0
  --conda-activate FILE      Default: /home/ps/anaconda3/bin/activate
  --metabolic-script FILE    METABOLIC-G.pl path
EOF
}

INPUT_DIR=""
OUTPUT_DIR=""
LOG_DIR=""
THREADS=100
MODULE_CUTOFF=0.75
KOFAM_DB=full
CONDA_ENV=METABOLIC_v4.0
CONDA_ACTIVATE=/home/ps/anaconda3/bin/activate
METABOLIC_SCRIPT=/home/ps/Research/data/working_dir/wangpudi/METABOLIC/METABOLIC-G.pl

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --log-dir) LOG_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --module-cutoff) MODULE_CUTOFF="$2"; shift 2 ;;
    --kofam-db) KOFAM_DB="$2"; shift 2 ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --conda-activate) CONDA_ACTIVATE="$2"; shift 2 ;;
    --metabolic-script) METABOLIC_SCRIPT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$INPUT_DIR" && -n "$OUTPUT_DIR" && -n "$LOG_DIR" ]] || {
  usage >&2
  exit 2
}
[[ -d "$INPUT_DIR" ]] || { echo "Input directory not found: $INPUT_DIR" >&2; exit 1; }
[[ -f "$CONDA_ACTIVATE" ]] || { echo "Conda activate script not found: $CONDA_ACTIVATE" >&2; exit 1; }
[[ -f "$METABOLIC_SCRIPT" ]] || { echo "METABOLIC-G.pl not found: $METABOLIC_SCRIPT" >&2; exit 1; }

INPUT_DIR=$(cd "$INPUT_DIR" && pwd -P)
OUTPUT_PARENT=$(dirname "$OUTPUT_DIR")
mkdir -p "$OUTPUT_PARENT" "$LOG_DIR"
OUTPUT_PARENT=$(cd "$OUTPUT_PARENT" && pwd -P)
OUTPUT_DIR="$OUTPUT_PARENT/$(basename "$OUTPUT_DIR")"
LOG_DIR=$(cd "$LOG_DIR" && pwd -P)

if [[ -e "$OUTPUT_DIR" ]] && find "$OUTPUT_DIR" -mindepth 1 -print -quit | grep -q .; then
  echo "Refusing to reuse a non-empty METABOLIC output directory: $OUTPUT_DIR" >&2
  exit 1
fi

mapfile -t FAA_FILES < <(
  find "$INPUT_DIR" -maxdepth 1 \( -type f -o -type l \) -name '*.faa' | sort
)
[[ ${#FAA_FILES[@]} -gt 0 ]] || { echo "No *.faa inputs found in $INPUT_DIR" >&2; exit 1; }
for faa in "${FAA_FILES[@]}"; do
  [[ -s "$faa" ]] || { echo "Empty or broken protein input: $faa" >&2; exit 1; }
done

# Some Conda activation hooks read optional variables before defining them.
# Keep strict mode for the workflow itself, but not while those hooks execute.
set +u
# shellcheck disable=SC1090
source "$CONDA_ACTIVATE" "$CONDA_ENV"
set -u

RUN_CONFIG="$LOG_DIR/run_config.tsv"
STDOUT_LOG="$LOG_DIR/METABOLIC_group_derep.stdout.log"
STARTED_AT=$(date -Iseconds)
{
  printf 'key\tvalue\n'
  printf 'started_at\t%s\n' "$STARTED_AT"
  printf 'input_dir\t%s\n' "$INPUT_DIR"
  printf 'input_MAG_count\t%s\n' "${#FAA_FILES[@]}"
  printf 'output_dir\t%s\n' "$OUTPUT_DIR"
  printf 'threads\t%s\n' "$THREADS"
  printf 'module_cutoff\t%s\n' "$MODULE_CUTOFF"
  printf 'kofam_db\t%s\n' "$KOFAM_DB"
  printf 'conda_env\t%s\n' "$CONDA_ENV"
  printf 'metabolic_script\t%s\n' "$METABOLIC_SCRIPT"
  printf 'metabolic_version\t4.0\n'
} > "$RUN_CONFIG"

COMMAND=(
  perl "$METABOLIC_SCRIPT"
  -t "$THREADS"
  -m-cutoff "$MODULE_CUTOFF"
  -in "$INPUT_DIR"
  -kofam-db "$KOFAM_DB"
  -o "$OUTPUT_DIR"
)
printf 'command=' | tee "$STDOUT_LOG"
printf ' %q' "${COMMAND[@]}" | tee -a "$STDOUT_LOG"
printf '\n' | tee -a "$STDOUT_LOG"

"${COMMAND[@]}" 2>&1 | tee -a "$STDOUT_LOG"

grep -q 'METABOLIC-G was done' "$OUTPUT_DIR/METABOLIC_log.log"
printf 'finished_at\t%s\n' "$(date -Iseconds)" >> "$RUN_CONFIG"
printf 'status\tcomplete\n' >> "$RUN_CONFIG"
echo "METABOLIC-G completed: $OUTPUT_DIR"
