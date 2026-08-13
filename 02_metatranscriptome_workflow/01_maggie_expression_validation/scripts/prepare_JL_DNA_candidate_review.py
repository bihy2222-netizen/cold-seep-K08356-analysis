#!/usr/bin/env python3
"""Build an auditable JL DNA candidate table without assigning IdrA function."""

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


def read_tbl(path):
    hits = {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith("#") and line.strip():
                fields = line.split()
                hits[fields[0]] = {"evalue": fields[4], "bitscore": fields[5]}
    return hits


def write_fasta(records, names, path):
    with open(path, "w") as handle:
        for name in names:
            if name in records:
                handle.write(f">{name}\n")
                seq = records[name]
                for start in range(0, len(seq), 80):
                    handle.write(seq[start : start + 80] + "\n")


def read_gff(path):
    genes = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = dict(
                item.split("=", 1) for item in fields[8].split(";") if "=" in item
            )
            local_id = attrs.get("ID", "")
            genes.append({
                "contig_id": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "local_id": local_id,
                "candidate_id": f"{fields[0]}_{local_id}" if local_id else "",
                "attributes": fields[8],
            })
    return genes


def write_neighborhoods(genes, candidate_names, path, flank=5):
    by_contig = {}
    for gene in genes:
        by_contig.setdefault(gene["contig_id"], []).append(gene)
    for contig_genes in by_contig.values():
        contig_genes.sort(key=lambda row: (row["start"], row["end"]))
    aliases = {}
    for contig_genes in by_contig.values():
        for index, gene in enumerate(contig_genes):
            aliases[gene["candidate_id"]] = (contig_genes, index)
            aliases[gene["local_id"]] = (contig_genes, index)
    with open(path, "w", newline="") as handle:
        fields = ["focal_candidate", "offset", "contig_id", "gene_id", "start", "end", "strand", "attributes"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for candidate in candidate_names:
            match = aliases.get(candidate)
            if match is None:
                continue
            contig_genes, focal_index = match
            left = max(0, focal_index - flank)
            right = min(len(contig_genes), focal_index + flank + 1)
            for index in range(left, right):
                gene = contig_genes[index]
                writer.writerow({
                    "focal_candidate": candidate,
                    "offset": index - focal_index,
                    "contig_id": gene["contig_id"],
                    "gene_id": gene["candidate_id"] or gene["local_id"],
                    "start": gene["start"],
                    "end": gene["end"],
                    "strand": gene["strand"],
                    "attributes": gene["attributes"],
                })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--combined", required=True)
    parser.add_argument("--combined-strict", required=True)
    parser.add_argument("--iria", required=True)
    parser.add_argument("--aioa", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--cds", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    combined = read_tbl(args.combined)
    strict = read_tbl(args.combined_strict)
    iria = read_tbl(args.iria)
    aioa = read_tbl(args.aioa)
    proteins = read_fasta(args.proteins)
    cds = read_fasta(args.cds)
    names = sorted(combined, key=lambda n: float(combined[n]["bitscore"]), reverse=True)

    fields = [
        "candidate_id", "combined_bitscore", "combined_evalue", "T640_status",
        "iriA_bitscore", "iriA_evalue", "aioA_bitscore", "aioA_evalue",
        "protein_length_aa", "tree_clade", "neighborhood_support", "final_class",
        "accept_for_personalized_reference", "review_notes",
    ]
    with open(outdir / "JL_0.1_K08356_candidate_review.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for name in names:
            writer.writerow({
                "candidate_id": name,
                "combined_bitscore": combined[name]["bitscore"],
                "combined_evalue": combined[name]["evalue"],
                "T640_status": "PASS" if name in strict else "exploratory",
                "iriA_bitscore": iria.get(name, {}).get("bitscore", ""),
                "iriA_evalue": iria.get(name, {}).get("evalue", ""),
                "aioA_bitscore": aioa.get(name, {}).get("bitscore", ""),
                "aioA_evalue": aioa.get(name, {}).get("evalue", ""),
                "protein_length_aa": len(proteins.get(name, "")),
                "tree_clade": "pending",
                "neighborhood_support": "pending",
                "final_class": "pending",
                "accept_for_personalized_reference": "no",
                "review_notes": "T640 is not an IdrA functional assignment",
            })
    write_fasta(proteins, names, outdir / "JL_0.1_K08356_candidates.faa")
    write_fasta(cds, names, outdir / "JL_0.1_K08356_candidates.CDS.fna")
    write_neighborhoods(
        read_gff(args.gff), names, outdir / "JL_0.1_K08356_candidate_neighborhood_coordinates.tsv"
    )


if __name__ == "__main__":
    main()
