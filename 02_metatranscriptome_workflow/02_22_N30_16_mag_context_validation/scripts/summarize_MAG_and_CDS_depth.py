#!/usr/bin/env python3
"""Stream all-bin depth by MAG and place IdrA among host-MAG CDS depths."""

import argparse
import bisect
import csv
import gzip
import statistics
from collections import Counter, defaultdict


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def parse_gff(path):
    by_contig = defaultdict(list)
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip(): continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS": continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            by_contig[fields[0]].append({"gene_id": attrs.get("ID", ""), "start": int(fields[3]),
                                         "end": int(fields[4]), "strand": fields[6]})
    return by_contig


def histogram_median(histogram, total):
    if not total: return 0
    targets = {(total - 1) // 2, total // 2}; observed = []
    cumulative = 0
    for depth, count in sorted(histogram.items()):
        previous = cumulative; cumulative += count
        for target in sorted(targets):
            if previous <= target < cumulative: observed.append(depth)
    return statistics.fmean(observed)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--depth", required=True); parser.add_argument("--gff", required=True)
    parser.add_argument("--bin-manifest", required=True); parser.add_argument("--strict-genes", required=True)
    parser.add_argument("--out-mag", required=True); parser.add_argument("--out-cds", required=True)
    parser.add_argument("--out-strict-cluster", required=True)
    parser.add_argument("--out-candidate-percentile", required=True)
    args = parser.parse_args()
    manifest = {row["reference_id"]: row for row in read_tsv(args.bin_manifest)}
    gff = parse_gff(args.gff)
    mag_stats = defaultdict(lambda: {"length":0, "sum":0, "ge1":0, "ge5":0, "ge10":0, "hist":Counter()})
    cds_rows, depths_by_mag = [], defaultdict(list)

    def finish_reference(ref, depths):
        if ref not in manifest: return
        meta = manifest[ref]; mag, contig = meta["MAG_id"], meta["original_contig_id"]
        stats = mag_stats[mag]; stats["length"] += len(depths); stats["sum"] += sum(depths)
        stats["ge1"] += sum(value >= 1 for value in depths)
        stats["ge5"] += sum(value >= 5 for value in depths)
        stats["ge10"] += sum(value >= 10 for value in depths); stats["hist"].update(depths)
        for gene in gff.get(contig, []):
            segment = depths[gene["start"]-1:gene["end"]]
            mean = statistics.fmean(segment) if segment else 0
            row = {"MAG_id": mag, "contig_id": contig, **gene, "CDS_length": len(segment),
                   "mean_depth": mean, "breadth_ge_1x": sum(x >= 1 for x in segment)/len(segment) if segment else 0}
            cds_rows.append(row); depths_by_mag[mag].append(mean)

    current_ref, values = None, []
    with gzip.open(args.depth, "rt") as handle:
        for line in handle:
            ref, _position, depth = line.rstrip().split("\t")[:3]
            if current_ref is not None and ref != current_ref:
                finish_reference(current_ref, values); values = []
            current_ref = ref; values.append(int(depth))
    if current_ref is not None: finish_reference(current_ref, values)

    mag_rows = []
    for mag, stats in sorted(mag_stats.items()):
        length = stats["length"]
        mag_rows.append({"MAG_id": mag, "MAG_length": length,
                         "breadth_ge_1x": stats["ge1"]/length if length else 0,
                         "breadth_ge_5x": stats["ge5"]/length if length else 0,
                         "breadth_ge_10x": stats["ge10"]/length if length else 0,
                         "mean_depth": stats["sum"]/length if length else 0,
                         "median_depth": histogram_median(stats["hist"], length)})
    write_tsv(args.out_mag, mag_rows, ["MAG_id","MAG_length","breadth_ge_1x","breadth_ge_5x","breadth_ge_10x","mean_depth","median_depth"])
    write_tsv(args.out_cds, cds_rows, ["MAG_id","contig_id","gene_id","start","end","strand","CDS_length","mean_depth","breadth_ge_1x"])

    cds_by_coordinates = {(row["MAG_id"], row["contig_id"], int(row["start"]), int(row["end"])): row for row in cds_rows}
    percentile_rows, strict_coverage_rows = [], []
    for strict in read_tsv(args.strict_genes):
        if not strict["MAG_id"]: continue
        key = (strict["MAG_id"], strict["contig_id"], int(strict["start"]), int(strict["end"]))
        observed = cds_by_coordinates.get(key)
        if not observed: continue
        strict_coverage_rows.append({"candidate_id": strict["candidate_id"], "role": strict["role"],
                                     "gene_id": strict["gene_id"], **observed})
        if strict["role"] != "IdrA": continue
        distribution = sorted(depths_by_mag[strict["MAG_id"]]); value = float(observed["mean_depth"])
        percentile = 100 * bisect.bisect_right(distribution, value) / len(distribution) if distribution else 0
        percentile_rows.append({"candidate_id": strict["candidate_id"], "MAG_id": strict["MAG_id"],
                                "IdrA_gene_id": strict["gene_id"], "IdrA_mean_depth": value,
                                "MAG_CDS_count": len(distribution), "IdrA_depth_percentile": percentile})
    write_tsv(args.out_candidate_percentile, percentile_rows,
              ["candidate_id","MAG_id","IdrA_gene_id","IdrA_mean_depth","MAG_CDS_count","IdrA_depth_percentile"])
    write_tsv(args.out_strict_cluster, strict_coverage_rows,
              ["candidate_id","role","gene_id","MAG_id","contig_id","start","end","strand","CDS_length","mean_depth","breadth_ge_1x"])


if __name__ == "__main__": main()
