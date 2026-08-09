#!/usr/bin/env python3
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(os.environ.get("ROOT", "/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra"))
REF_DIR = ROOT / "05_idra_reference"


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path, errors="replace") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


def write_fasta(path, records):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def normalize_level(value):
    value = (value or "").lower()
    if "strict" in value or "dirm" in value:
        return "strict"
    if "extended" in value or "partial" in value or "idr" in value:
        return "extended"
    if "competitive" in value or "aio" in value or "dmsor" in value or "arr" in value or "nar" in value:
        return "competitive_only"
    return "exclude"


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 02_prepare_idra_cds_references.py 05_idra_reference/idra_reference_metadata.tsv")
    metadata = Path(sys.argv[1])
    if not metadata.exists():
        raise SystemExit(f"Missing reference metadata: {metadata}")
    REF_DIR.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(metadata), delimiter="\t"))
    strict, extended, competitive = [], [], []
    seen = set()
    qc_rows = []

    for row in rows:
        gene_id = row.get("gene_id") or row.get("id") or row.get("sequence_id")
        cds_path = Path(row.get("original_CDS_path") or row.get("CDS_path") or row.get("cds_path") or "")
        level = normalize_level(row.get("reference_level") or row.get("clade") or row.get("neighbor_category"))
        if not gene_id or not cds_path.exists():
            row["reference_qc_status"] = "missing_gene_id_or_CDS_path"
            qc_rows.append(row)
            continue
        seqs = read_fasta(cds_path)
        seq = seqs.get(gene_id)
        if seq is None and len(seqs) == 1:
            seq = next(iter(seqs.values()))
        if seq is None:
            row["reference_qc_status"] = "gene_id_not_found_in_CDS_fasta"
            qc_rows.append(row)
            continue
        seq = re.sub(r"\s+", "", seq).upper().replace("U", "T")
        if gene_id in seen:
            row["reference_qc_status"] = "duplicate_gene_id_skipped"
            qc_rows.append(row)
            continue
        seen.add(gene_id)
        header = "|".join([
            gene_id,
            row.get("MAG_id", ""),
            row.get("clade", ""),
            row.get("neighbor_category", ""),
            row.get("host_taxonomy", ""),
        ]).replace(" ", "_")
        rec = (header, seq)
        if level == "strict":
            strict.append(rec)
            extended.append(rec)
            competitive.append(rec)
        elif level == "extended":
            extended.append(rec)
            competitive.append(rec)
        elif level == "competitive_only":
            competitive.append(rec)
        row["reference_qc_status"] = "included_" + level
        row["CDS_length_observed"] = str(len(seq))
        qc_rows.append(row)

    write_fasta(REF_DIR / "idra_strict_20_CDS.fna", strict)
    write_fasta(REF_DIR / "idra_extended_CDS.fna", extended)
    write_fasta(REF_DIR / "DMSOR_competitive_CDS.fna", competitive)

    fields = list(dict.fromkeys(k for row in qc_rows for k in row.keys()))
    with open(REF_DIR / "idra_reference_metadata.qc.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(qc_rows)

    for fasta, index in [
        ("idra_strict_20_CDS.fna", "idra_strict_20"),
        ("idra_extended_CDS.fna", "idra_extended"),
        ("DMSOR_competitive_CDS.fna", "DMSOR_competitive"),
    ]:
        fasta_path = REF_DIR / fasta
        if fasta_path.stat().st_size == 0:
            print(f"WARNING empty reference: {fasta_path}")
            continue
        subprocess.run(["bowtie2-build", str(fasta_path), str(REF_DIR / index)], check=True)
        subprocess.run(["samtools", "faidx", str(fasta_path)], check=True)

    print(f"strict_CDS={len(strict)}")
    print(f"extended_CDS={len(extended)}")
    print(f"competitive_CDS={len(competitive)}")
    print(f"reference_dir={REF_DIR}")


if __name__ == "__main__":
    main()
