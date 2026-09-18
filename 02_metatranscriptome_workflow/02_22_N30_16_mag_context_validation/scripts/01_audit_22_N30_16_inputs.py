#!/usr/bin/env python3
"""Inventory existing 22_N30_16 inputs without modifying source results."""

import argparse
import csv
import hashlib
from pathlib import Path


PATTERNS = {
    "raw_DNA_R1": ["**/SRR35213524_1.fastq.gz"],
    "raw_DNA_R2": ["**/SRR35213524_2.fastq.gz"],
    "clean_DNA_R1": ["**/22_N30_16*R1*fastp*f*q.gz", "**/SRR35213524*clean*R1*f*q.gz"],
    "clean_DNA_R2": ["**/22_N30_16*R2*fastp*f*q.gz", "**/SRR35213524*clean*R2*f*q.gz"],
    "assembly": ["22_N30_16/**/final.contigs.fa"],
    "prodigal_FAA": ["22_N30_16/**/*protein*.faa", "22_N30_16/**/*.proteins.faa"],
    "prodigal_FNA": ["22_N30_16/**/*gene*.fna", "22_N30_16/**/*.CDS.fna"],
    "prodigal_GFF": ["22_N30_16/**/*.gff"],
    "combined_HMM_T640": ["22_N30_16/**/*combined*T640*.tbl*"],
    "combined_HMM_all": ["22_N30_16/**/*combined*all*.tbl*"],
    "iriA_HMM": ["22_N30_16/**/*iriA*.tbl*"],
    "aioA_HMM": ["22_N30_16/**/*aioA*.tbl*"],
    "treefile": ["**/*22_N30_16*.treefile", "**/*candidate*.treefile"],
    "label_map": ["**/*label*map*.tsv"],
    "bin_FASTA": ["22_N30_16/**/*bin*.fa", "22_N30_16/**/*bin*.fasta"],
    "contig_to_bin": ["22_N30_16/**/*contig*bin*.tsv"],
    "neighborhood": ["22_N30_16/**/*neigh*.tsv"],
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    checksum_lines = []
    for file_class, patterns in PATTERNS.items():
        matches = set()
        for pattern in patterns:
            matches.update(path for path in args.pilot_root.glob(pattern) if path.is_file())
        status = "MISSING" if not matches else "UNIQUE" if len(matches) == 1 else "AMBIGUOUS"
        if not matches:
            rows.append({"file_class": file_class, "status": status, "path": "", "bytes": "", "sha256": ""})
        for path in sorted(matches):
            digest = sha256(path)
            rows.append({"file_class": file_class, "status": status, "path": str(path.resolve()),
                         "bytes": path.stat().st_size, "sha256": digest})
            checksum_lines.append(f"{digest}  {path.resolve()}")
    with (args.outdir / "input_file_inventory.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    (args.outdir / "checksums.sha256").write_text("\n".join(checksum_lines) + "\n")
    unresolved = sorted({r["file_class"] for r in rows if r["status"] != "UNIQUE"})
    (args.outdir / "AUDIT_STATUS.txt").write_text(
        "Manual path resolution required for:\n" + "\n".join(unresolved) + "\n"
        if unresolved else "All audited file classes have one match.\n"
    )


if __name__ == "__main__":
    main()

