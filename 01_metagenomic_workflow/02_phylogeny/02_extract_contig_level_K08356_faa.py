#!/usr/bin/env python3
"""Extract contig-level K08356/AioA protein sequences from KOfamScan best hits.

Original server context:
  /home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/CDS

This script reproduces the original contig-level extraction used to create:
  aio-A_contig_level/*.K08356.faa
  aio-A_contig_level/K08356_geneid_contigid_summary.tsv

It reads per-sample *.cdhit.kofamscan.best.txt files, keeps the requested KO
from the best/reliable KOfamScan hits, and extracts matching proteins from the
corresponding *.cdhit.cds.faa files.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def get_contig_id(gene_id: str) -> str:
    """Convert a Prodigal gene ID to its parent contig ID."""
    if "_" in gene_id:
        return gene_id.rsplit("_", 1)[0]
    return gene_id


def parse_kofam_best_file(kofam_file: Path, target_ko: str) -> list[dict[str, str]]:
    """Parse one *.cdhit.kofamscan.best.txt file for the target KO."""
    records: list[dict[str, str]] = []

    with kofam_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            if "\t" in line:
                parts = line.split("\t")
            else:
                parts = re.split(r"\s+", line, maxsplit=6)

            if len(parts) < 3:
                continue

            if parts[0] == "*":
                gene_id = parts[1]
                ko_id = parts[2]
                threshold = parts[3] if len(parts) > 3 else ""
                score = parts[4] if len(parts) > 4 else ""
                evalue = parts[5] if len(parts) > 5 else ""
                annotation = parts[6] if len(parts) > 6 else ""
            else:
                gene_id = parts[0]
                ko_id = parts[1]
                threshold = parts[2] if len(parts) > 2 else ""
                score = parts[3] if len(parts) > 3 else ""
                evalue = parts[4] if len(parts) > 4 else ""
                annotation = parts[5] if len(parts) > 5 else ""

            if ko_id == target_ko:
                records.append(
                    {
                        "gene_id": gene_id,
                        "contig_id": get_contig_id(gene_id),
                        "ko_id": ko_id,
                        "kofam_threshold": threshold,
                        "kofam_score": score,
                        "evalue": evalue,
                        "annotation": annotation.strip('"'),
                    }
                )

    return records


def read_fasta(fasta_file: Path) -> dict[str, str]:
    """Read a FASTA file as {sequence_id: sequence}."""
    seqs: dict[str, str] = {}
    current_id: str | None = None
    seq_lines: list[str] = []

    with fasta_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(seq_lines)
                current_id = line[1:].strip().split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())

        if current_id is not None:
            seqs[current_id] = "".join(seq_lines)

    return seqs


def write_fasta(
    records: list[dict[str, str]],
    seqs: dict[str, str],
    sample_name: str,
    out_faa: Path,
    missing_records: list[dict[str, str]],
) -> int:
    """Write extracted proteins with sample-prefixed sequence IDs."""
    n_written = 0

    with out_faa.open("w", encoding="utf-8") as out:
        for rec in records:
            gene_id = rec["gene_id"]
            if gene_id not in seqs:
                missing_records.append(
                    {
                        "sample": sample_name,
                        "gene_id": gene_id,
                        "contig_id": rec["contig_id"],
                        "ko_id": rec["ko_id"],
                    }
                )
                continue

            out.write(f">{sample_name}-{gene_id}\n")
            seq = seqs[gene_id]
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")
            n_written += 1

    if n_written == 0 and out_faa.exists() and out_faa.stat().st_size == 0:
        out_faa.unlink()

    return n_written


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract proteins for one KO from *.cdhit.kofamscan.best.txt and "
            "matching *.cdhit.cds.faa files."
        )
    )
    parser.add_argument("-k", "--ko", default="K08356", help="Target KO ID.")
    parser.add_argument(
        "-i",
        "--input-dir",
        default=".",
        help="Directory containing *.cdhit.kofamscan.best.txt and *.cdhit.cds.faa.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="aio-A_contig_level",
        help="Directory for per-sample KO FASTA files and summary tables.",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = input_dir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    target_ko = args.ko
    kofam_files = sorted(input_dir.glob("*.cdhit.kofamscan.best.txt"))

    if not kofam_files:
        raise SystemExit(
            f"ERROR: no *.cdhit.kofamscan.best.txt files found in {input_dir}"
        )

    summary_file = output_dir / f"{target_ko}_geneid_contigid_summary.tsv"
    missing_file = output_dir / f"{target_ko}_missing_geneid_in_faa.tsv"

    all_summary_records: list[dict[str, str]] = []
    missing_records: list[dict[str, str]] = []

    print(f"Target KO: {target_ko}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Found {len(kofam_files)} KOfamScan best files.\n")

    for kofam_file in kofam_files:
        sample_name = kofam_file.name.replace(".cdhit.kofamscan.best.txt", "")
        faa_file = input_dir / f"{sample_name}.cdhit.cds.faa"
        out_faa = output_dir / f"{sample_name}.{target_ko}.faa"

        records = parse_kofam_best_file(kofam_file, target_ko)
        if not records:
            print(f"[SKIP] {sample_name}: no {target_ko} best hits.")
            continue

        for rec in records:
            tmp = rec.copy()
            tmp["sample"] = sample_name
            tmp["kofam_file"] = kofam_file.name
            tmp["faa_file"] = faa_file.name
            all_summary_records.append(tmp)

        if not faa_file.exists():
            print(
                f"[WARN] {sample_name}: {len(records)} hits, but {faa_file.name} "
                "was not found."
            )
            for rec in records:
                missing_records.append(
                    {
                        "sample": sample_name,
                        "gene_id": rec["gene_id"],
                        "contig_id": rec["contig_id"],
                        "ko_id": rec["ko_id"],
                    }
                )
            continue

        seqs = read_fasta(faa_file)
        n_written = write_fasta(records, seqs, sample_name, out_faa, missing_records)
        print(f"[OK] {sample_name}: {n_written}/{len(records)} sequences -> {out_faa}")

    header = [
        "sample",
        "gene_id",
        "contig_id",
        "ko_id",
        "kofam_threshold",
        "kofam_score",
        "evalue",
        "annotation",
        "kofam_file",
        "faa_file",
    ]
    with summary_file.open("w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for rec in all_summary_records:
            out.write("\t".join(str(rec.get(col, "")) for col in header) + "\n")

    if missing_records:
        with missing_file.open("w", encoding="utf-8") as out:
            missing_header = ["sample", "gene_id", "contig_id", "ko_id"]
            out.write("\t".join(missing_header) + "\n")
            for rec in missing_records:
                out.write("\t".join(str(rec.get(col, "")) for col in missing_header) + "\n")
        print(f"\nMissing gene IDs written to: {missing_file}")

    print(f"\nSummary table written to: {summary_file}")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
