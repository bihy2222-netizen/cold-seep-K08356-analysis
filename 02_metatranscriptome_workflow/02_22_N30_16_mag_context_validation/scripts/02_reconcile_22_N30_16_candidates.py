#!/usr/bin/env python3
"""Reconcile exact Prodigal/HMM IDs and extract candidate neighborhoods."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def fasta(path):
    records, descriptions, name, seq = {}, {}, None, []
    for line in Path(path).open(errors="replace"):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(seq)
            description = line[1:]
            name, seq = description.split()[0], []
            descriptions[name] = description
        elif name is not None:
            seq.append(line)
    if name is not None:
        records[name] = "".join(seq)
    return records, descriptions


def tbl(path):
    hits = {}
    for line in Path(path).open(errors="replace"):
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        hits[fields[0]] = {"evalue": fields[4], "bitscore": float(fields[5])}
    return hits


def domtbl(path):
    hits = defaultdict(list)
    for line in Path(path).open(errors="replace"):
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        target_length = int(fields[2])
        alignment_start, alignment_end = int(fields[17]), int(fields[18])
        coverage = max(0, alignment_end - alignment_start + 1) / target_length if target_length else 0
        hits[fields[0]].append((coverage, alignment_start, alignment_end))
    return hits


def parse_attributes(text):
    return dict(item.split("=", 1) for item in text.split(";") if "=" in item)


def gff(path):
    rows, by_contig = [], defaultdict(list)
    for line in Path(path).open(errors="replace"):
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip().split("\t")
        if len(fields) != 9 or fields[2] != "CDS":
            continue
        attrs = parse_attributes(fields[8])
        row = {
            "gff_id": attrs.get("ID", ""), "contig_id": fields[0],
            "start": int(fields[3]), "end": int(fields[4]), "strand": fields[6],
            "attributes": fields[8],
        }
        rows.append(row)
        by_contig[fields[0]].append(row)
    for genes in by_contig.values():
        genes.sort(key=lambda item: (item["start"], item["end"]))
    return rows, by_contig


def prodigal_coordinates(descriptions):
    """Read Prodigal '# start # end # strand' headers without fuzzy IDs."""
    coords = {}
    for sequence_id, description in descriptions.items():
        pieces = [piece.strip() for piece in description.split("#")]
        if len(pieces) < 4:
            continue
        try:
            start, end, strand_code = int(pieces[1]), int(pieces[2]), int(pieces[3])
        except ValueError:
            continue
        contig_id = sequence_id.rsplit("_", 1)[0]
        strand = "+" if strand_code == 1 else "-"
        coords[(contig_id, start, end, strand)] = sequence_id
    return coords


def reconcile_gff_ids(gff_rows, protein_descriptions):
    by_coord = prodigal_coordinates(protein_descriptions)
    crosswalk = []
    for row in gff_rows:
        key = (row["contig_id"], row["start"], row["end"], row["strand"])
        protein_id = by_coord.get(key)
        method = "exact_coordinate_and_strand" if protein_id else "unresolved"
        if not protein_id and row["gff_id"] in protein_descriptions:
            protein_id, method = row["gff_id"], "exact_identifier"
        crosswalk.append({**row, "protein_id": protein_id or "", "crosswalk_method": method})
    return crosswalk


def write_fasta(records, names, path):
    with Path(path).open("w") as handle:
        for name in names:
            if name not in records:
                continue
            handle.write(f">{name}\n")
            for start in range(0, len(records[name]), 80):
                handle.write(records[name][start:start + 80] + "\n")


def main():
    parser = argparse.ArgumentParser()
    for key in ("combined", "combined_strict", "combined_domtbl", "iria", "aioa",
                "faa", "fna", "gff", "contigs", "outdir"):
        parser.add_argument("--" + key.replace("_", "-"), required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    combined, strict = tbl(args.combined), tbl(args.combined_strict)
    iria, aioa, domains = tbl(args.iria), tbl(args.aioa), domtbl(args.combined_domtbl)
    proteins, protein_descriptions = fasta(args.faa)
    cds, _ = fasta(args.fna)
    contigs, _ = fasta(args.contigs)
    gff_rows, _ = gff(args.gff)
    crosswalk = reconcile_gff_ids(gff_rows, protein_descriptions)
    by_protein = {row["protein_id"]: row for row in crosswalk if row["protein_id"]}
    by_contig = defaultdict(list)
    for row in crosswalk:
        by_contig[row["contig_id"]].append(row)
    for genes in by_contig.values():
        genes.sort(key=lambda item: (item["start"], item["end"]))

    names = sorted(strict, key=lambda item: strict[item]["bitscore"], reverse=True)
    if len(names) != 10:
        raise SystemExit(f"Expected 10 T640 candidates for 22_N30_16; found {len(names)}")
    rows, neighborhoods = [], []
    for name in names:
        if name not in proteins or name not in cds or name not in by_protein:
            raise SystemExit(f"Exact/coordinate ID reconciliation failed for {name}; no fuzzy join was attempted")
        gene = by_protein[name]
        contig_length = len(contigs.get(gene["contig_id"], ""))
        edge_distance = min(gene["start"] - 1, contig_length - gene["end"]) if contig_length else ""
        best_domain = max(domains.get(name, [(0, "", "")]), key=lambda item: item[0])
        rows.append({
            "sample": "22_N30_16", "contig_id": gene["contig_id"], "protein_id": name,
            "CDS_id": name, "strand": gene["strand"], "CDS_start": gene["start"],
            "CDS_end": gene["end"], "CDS_length": len(cds[name]),
            "protein_length": len(proteins[name].rstrip("*")),
            "combined_bitscore": combined.get(name, {}).get("bitscore", ""),
            "iriA_bitscore": iria.get(name, {}).get("bitscore", ""),
            "aioA_bitscore": aioa.get(name, {}).get("bitscore", ""),
            "HMM_alignment_coverage": round(best_domain[0], 6),
            "HMM_alignment_start_aa": best_domain[1], "HMM_alignment_end_aa": best_domain[2],
            "contig_length": contig_length, "contig_edge_distance": edge_distance,
            "tree_position": "pending", "nearest_reference": "pending", "tree_support": "pending",
            "final_phylogenetic_class": "pending",
        })
        genes = by_contig[gene["contig_id"]]
        index = next(i for i, item in enumerate(genes) if item["protein_id"] == name)
        edge_censored = index < 10 or len(genes) - index - 1 < 10
        for neighbor_index in range(max(0, index - 10), min(len(genes), index + 11)):
            neighbor = genes[neighbor_index]
            neighborhoods.append({
                "candidate_id": name, "relative_ORF": neighbor_index - index,
                "contig_id": neighbor["contig_id"], "neighbor_id": neighbor["protein_id"],
                "gff_id": neighbor["gff_id"], "crosswalk_method": neighbor["crosswalk_method"],
                "start": neighbor["start"], "end": neighbor["end"], "strand": neighbor["strand"],
                "contig_length": contig_length, "edge_censored": "yes" if edge_censored else "no",
                "neighbor_family": "pending", "neighbor_annotation_accepted": "no",
                "neighbor_bitscore": "", "neighbor_reference": "",
            })

    for path, data in [
        (outdir / "candidate_id_crosswalk.tsv", crosswalk),
        (outdir / "22_N30_16_candidate_review.tsv", rows),
        (outdir / "22_N30_16_candidate_neighborhood.tsv", neighborhoods),
    ]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(data)
    write_fasta(proteins, names, outdir / "22_N30_16_T640_candidates.faa")
    write_fasta(cds, names, outdir / "22_N30_16_T640_candidates.CDS.fna")
    write_fasta(proteins, sorted({row["neighbor_id"] for row in neighborhoods if row["neighbor_id"]}),
                outdir / "22_N30_16_candidate_neighborhood_proteins.faa")
    write_fasta(contigs, sorted({row["contig_id"] for row in rows}), outdir / "22_N30_16_candidate_contigs.fna")


if __name__ == "__main__":
    main()
