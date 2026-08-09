#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra}"
SEARCH_ROOT="${SEARCH_ROOT:-/home/ps/ps1/data/bihongyu/cold_seep}"

mkdir -p "$ROOT"/{00_metadata,01_sra,02_raw_fastq,03_fastp,04_non_rRNA,05_idra_reference,06_bowtie2,07_assembly,08_prodigal,09_hmmsearch,10_expression,11_summary,logs,tmp,scripts}
cd "$ROOT"

{
  echo "# Resource check"
  date
  pwd
  df -h
  df -i
  free -h || true
  nproc || true
} > 00_metadata/resource_check.txt

TOOLS=(
  esearch efetch prefetch fasterq-dump fastp sortmerna bowtie2 bowtie2-build
  samtools megahit prodigal hmmsearch hmmpress seqkit coverm salmon pigz
  md5sum gzip awk sed find python3
)

: > 00_metadata/software_versions.txt
for tool in "${TOOLS[@]}"; do
  {
    echo "## $tool"
    if command -v "$tool" >/dev/null 2>&1; then
      command -v "$tool"
      "$tool" --version 2>&1 | head -n 5 || true
    else
      echo "MISSING"
    fi
    echo
  } >> 00_metadata/software_versions.txt
done

find "$SEARCH_ROOT" \
  -type f \
  \( -iname 'idra.hmm' -o -iname '*idra*.hmm' -o -iname '*iriA*.hmm' -o -iname '*aioA*.hmm' \) \
  2>/dev/null \
  | sort > 00_metadata/hmm_candidate_paths.txt

: > 00_metadata/hmm_candidate_summary.tsv
printf "path\tsize_bytes\tmtime\tsha256\tNAME\tACC\tLENG\tNSEQ\tmodel_note\n" > 00_metadata/hmm_candidate_summary.tsv
while IFS= read -r hmm; do
  [ -s "$hmm" ] || continue
  name="$(awk '$1=="NAME"{print $2; exit}' "$hmm" 2>/dev/null || true)"
  acc="$(awk '$1=="ACC"{print $2; exit}' "$hmm" 2>/dev/null || true)"
  leng="$(awk '$1=="LENG"{print $2; exit}' "$hmm" 2>/dev/null || true)"
  nseq="$(awk '$1=="NSEQ"{print $2; exit}' "$hmm" 2>/dev/null || true)"
  sha="$(sha256sum "$hmm" | awk '{print $1}')"
  size="$(stat -c '%s' "$hmm")"
  mtime="$(stat -c '%y' "$hmm")"
  note="candidate"
  case "$(basename "$hmm" | tr '[:upper:]' '[:lower:]')" in
    idra.hmm) note="preferred_name_idra_hmm" ;;
    *combined*|*aioa*) note="combined_or_aioA_related_not_standalone_idrA" ;;
    *iria*) note="iriA_idrA_related_candidate" ;;
  esac
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$hmm" "$size" "$mtime" "$sha" "$name" "$acc" "$leng" "$nseq" "$note" >> 00_metadata/hmm_candidate_summary.tsv
done < 00_metadata/hmm_candidate_paths.txt

find "$SEARCH_ROOT" \
  -type f \
  \( -iname '*K08356*.faa' -o -iname '*K08356*.fna' -o -iname '*idra*.faa' -o -iname '*idra*.fna' -o -iname '*candidate*.faa' -o -iname '*candidate*.fna' \) \
  2>/dev/null \
  | sort > 00_metadata/idra_reference_candidate_paths.txt

cp 00_metadata/software_versions.txt 11_summary/software_versions.txt
echo "Wrote 00_metadata/resource_check.txt, software_versions.txt, hmm_candidate_summary.tsv, idra_reference_candidate_paths.txt"
