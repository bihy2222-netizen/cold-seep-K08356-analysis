#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${1:?Usage: bash 02_download_fastq_from_manifest.sh MANIFEST.tsv OUT_DIR}"
OUT_DIR="${2:?Usage: bash 02_download_fastq_from_manifest.sh MANIFEST.tsv OUT_DIR}"
ROOT="${ROOT:-/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra}"

if [ "${CONFIRM_QDN_RNA_DOWNLOAD:-}" != "YES" ]; then
  echo "Refusing to download. Set CONFIRM_QDN_RNA_DOWNLOAD=YES after checking the manifest." >&2
  exit 2
fi

mkdir -p "$OUT_DIR" "$ROOT/logs" "$ROOT/00_metadata"
STATUS="$ROOT/00_metadata/download_status.tsv"
MD5_OUT="$ROOT/00_metadata/fastq_md5.tsv"

if [ ! -s "$STATUS" ]; then
  printf "timestamp\trun\tsample\tfile\tstatus\tnote\n" > "$STATUS"
fi
if [ ! -s "$MD5_OUT" ]; then
  printf "run\tsample\tfile\tmd5_expected\tmd5_observed\tgzip_test\n" > "$MD5_OUT"
fi

python3 - "$MANIFEST" "$OUT_DIR" "$STATUS" "$MD5_OUT" <<'PY'
import csv
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

manifest, out_dir, status_path, md5_path = map(Path, sys.argv[1:5])
rows = list(csv.DictReader(open(manifest), delimiter="\t"))

def parts(value):
    return [x for x in str(value or "").split(";") if x]

def log_status(run, sample, file_name, status, note=""):
    with open(status_path, "a") as h:
        h.write("\t".join([datetime.now().isoformat(timespec="seconds"), run, sample, file_name, status, note]) + "\n")

def md5sum(path):
    out = subprocess.check_output(["md5sum", str(path)], text=True)
    return out.split()[0]

for row in rows:
    run = row.get("Run") or row.get("run_accession") or row.get("ENA_run_accession")
    sample = row.get("SampleName") or row.get("Sample_Name") or row.get("ENA_sample_accession") or "sample"
    urls = parts(row.get("FASTQ_URL") or row.get("ENA_fastq_ftp"))
    urls = ["https://" + u if u.startswith("ftp.sra.ebi.ac.uk") else u for u in urls]
    md5s = parts(row.get("FASTQ_MD5") or row.get("ENA_fastq_md5"))
    sample_dir = out_dir / f"{sample}_{run}"
    sample_dir.mkdir(parents=True, exist_ok=True)
    for idx, url in enumerate(urls, start=1):
        expected = md5s[idx - 1] if idx - 1 < len(md5s) else ""
        fname = url.rsplit("/", 1)[-1]
        final = sample_dir / fname
        part = sample_dir / (fname + ".part")
        if final.exists() and final.stat().st_size > 0:
            observed = md5sum(final)
            if expected and observed == expected:
                subprocess.run(["gzip", "-t", str(final)], check=True)
                log_status(run, sample, fname, "existing_md5_gzip_ok")
                with open(md5_path, "a") as h:
                    h.write("\t".join([run, sample, str(final), expected, observed, "PASS"]) + "\n")
                continue
            log_status(run, sample, fname, "existing_md5_failed", observed)
        if part.exists():
            part.unlink()
        log_status(run, sample, fname, "download_started", url)
        with open(sample_dir / f"{fname}.wget.log", "w") as log:
            subprocess.run(["wget", "-O", str(part), url], stdout=log, stderr=subprocess.STDOUT, check=True)
        observed = md5sum(part)
        if expected and observed != expected:
            log_status(run, sample, fname, "md5_failed", observed)
            raise SystemExit(f"MD5 failed for {part}: expected {expected}, observed {observed}")
        subprocess.run(["gzip", "-t", str(part)], check=True)
        part.rename(final)
        log_status(run, sample, fname, "download_md5_gzip_ok")
        with open(md5_path, "a") as h:
            h.write("\t".join([run, sample, str(final), expected, observed, "PASS"]) + "\n")
PY

echo "Download workflow finished: $(date)"
