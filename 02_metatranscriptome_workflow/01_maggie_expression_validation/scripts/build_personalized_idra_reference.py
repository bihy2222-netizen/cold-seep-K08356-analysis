#!/usr/bin/env python3
"""Extract only tree- and synteny-approved sample-specific IdrA CDS."""

import argparse
import csv
from pathlib import Path


def read_fasta(path):
    records = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks)
    return records


def write_fasta(records, names, path):
    with open(path, "w") as handle:
        for name in names:
            handle.write(f">JL_0.1_personalized|{name}\n")
            seq = records[name]
            for start in range(0, len(seq), 80):
                handle.write(seq[start : start + 80] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--review-tsv", required=True)
    parser.add_argument("--cds", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cds = read_fasta(args.cds)
    proteins = read_fasta(args.proteins)

    accepted = []
    with open(args.review_tsv, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row["final_class"] == "strict_IdrA_associated"
                    and row["neighborhood_support"] == "complete_DIRM_like"
                    and row["accept_for_personalized_reference"].lower() == "yes"):
                accepted.append(row["candidate_id"])
    missing = [name for name in accepted if name not in cds or name not in proteins]
    if missing:
        raise SystemExit("Accepted IDs missing from CDS/protein FASTA: " + ", ".join(missing))
    if not accepted:
        raise SystemExit("No tree- and synteny-approved strict IdrA candidates were accepted")

    write_fasta(cds, accepted, outdir / "JL_0.1_personalized_strict_IdrA_CDS.fna")
    write_fasta(proteins, accepted, outdir / "JL_0.1_personalized_strict_IdrA_proteins.faa")
    with open(outdir / "JL_0.1_personalized_reference_manifest.tsv", "w") as handle:
        handle.write("candidate_id\tcds_length_nt\tprotein_length_aa\tstatus\n")
        for name in accepted:
            handle.write(f"{name}\t{len(cds[name])}\t{len(proteins[name])}\taccepted\n")


if __name__ == "__main__":
    main()

