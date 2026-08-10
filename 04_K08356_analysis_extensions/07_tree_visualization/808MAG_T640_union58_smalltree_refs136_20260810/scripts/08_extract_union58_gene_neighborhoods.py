#!/usr/bin/env python3
"""Extract +/-10 CDS neighborhoods for the 58 MAG candidates from per-MAG Prodigal proteins."""

import argparse
import csv
import re
from pathlib import Path


HEADER_RE = re.compile(r"^>(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#\s+(.*)$")
ID_RE = re.compile(r"^(.+)-k141_(\d+)_(\d+)$")


def fasta_records(path):
    header = None
    seq = []
    for line in path.open(errors="replace"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header, seq = line, []
        elif header is not None:
            seq.append(line.strip().rstrip("*"))
    if header is not None:
        yield header, "".join(seq)


def parse_header(header):
    match = HEADER_RE.match(header)
    if not match:
        return None
    seq_id, start, end, strand, attrs = match.groups()
    ident = ID_RE.match(seq_id)
    if not ident:
        return None
    contig, contig_number, gene_index = ident.groups()
    return {
        "gene_id": seq_id,
        "contig_id": f"{contig}-k141_{contig_number}",
        "gene_index": int(gene_index),
        "start": int(start),
        "end": int(end),
        "strand": "+" if strand == "1" else "-",
        "attrs": attrs,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-table", type=Path, required=True)
    parser.add_argument("--protein-dir", type=Path, required=True)
    parser.add_argument("--out-long", type=Path, required=True)
    parser.add_argument("--out-faa", type=Path, required=True)
    args = parser.parse_args()

    with args.candidate_table.open(newline="") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))

    long_rows = []
    query_records = {}
    for candidate in candidates:
        target_id = candidate["protein_ID"]
        mag_id = candidate["MAG_ID"]
        target_match = ID_RE.match(target_id)
        if not target_match:
            raise SystemExit(f"Cannot parse target ID: {target_id}")
        contig_prefix, contig_number, target_index = target_match.groups()
        target_contig = f"{contig_prefix}-k141_{contig_number}"
        protein_path = args.protein_dir / f"{mag_id}.faa"
        if not protein_path.exists():
            raise SystemExit(f"Missing per-MAG protein FASTA: {protein_path}")

        found_target = False
        for header, sequence in fasta_records(protein_path):
            parsed = parse_header(header)
            if not parsed or parsed["contig_id"] != target_contig:
                continue
            relative = parsed["gene_index"] - int(target_index)
            if abs(relative) > 10:
                continue
            found_target |= relative == 0
            query_records[parsed["gene_id"]] = sequence
            long_rows.append({
                "target_id": target_id,
                "MAG_ID": mag_id,
                "contig_id": target_contig,
                "target_gene_index": target_index,
                "neighbor_gene_index": parsed["gene_index"],
                "relative_position": relative,
                "role": "target" if relative == 0 else ("upstream" if relative < 0 else "downstream"),
                "neighbor_raw_id": parsed["gene_id"],
                "start": parsed["start"],
                "end": parsed["end"],
                "strand": parsed["strand"],
                "prodigal_attrs": parsed["attrs"],
            })
        if not found_target:
            raise SystemExit(f"Target not found in per-MAG protein FASTA: {target_id}")

    long_rows.sort(key=lambda row: (row["target_id"], int(row["neighbor_gene_index"])))
    args.out_long.parent.mkdir(parents=True, exist_ok=True)
    fields = list(long_rows[0])
    with args.out_long.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(long_rows)
    with args.out_faa.open("w") as handle:
        for gene_id, sequence in query_records.items():
            handle.write(f">{gene_id}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start:start + 70] + "\n")

    print(f"candidates={len(candidates)}")
    print(f"neighborhood_rows={len(long_rows)}")
    print(f"query_proteins={len(query_records)}")


if __name__ == "__main__":
    main()
