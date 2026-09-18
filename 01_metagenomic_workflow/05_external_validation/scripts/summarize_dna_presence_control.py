#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path


EXCLUDE_PRIMARY = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
DEPTH_PROFILES = {
    "raw_mapped": ["-g", "0x700", "-G", "0x4", "-Q", "0"],
    "primary_deduplicated": ["-G", "0x800", "-Q", "0"],
    "primary_MAPQ20_deduplicated": ["-G", "0x800", "-Q", "20"],
    "primary_MAPQ30_deduplicated": ["-G", "0x800", "-Q", "30"],
}


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames=None):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def make_regions(classification, neighborhoods, sample):
    candidates = {row["candidate_id"]: row for row in classification if row["sample"] == sample}
    grouped = defaultdict(list)
    for row in neighborhoods:
        if row["sample"] == sample:
            grouped[row["candidate_id"]].append(row)
    regions, clusters = [], {}
    for candidate_id, candidate in candidates.items():
        focal = next(row for row in grouped[candidate_id] if int(row["relative_ORF"]) == 0)
        base = {
            "sample": sample, "candidate_id": candidate_id, "candidate_protein_id": candidate["protein_id"],
            "contig_id": focal["contig_id"], "final_integrated_class": candidate["final_integrated_class"],
            "bin_id": candidate["bin_id"], "GTDB_family": candidate["GTDB_family"],
        }
        focal_region = {**base, "gene_role": "IdrA_candidate", "gene_id": focal["neighbor_id"],
                        "start": int(focal["start"]), "end": int(focal["end"]), "strand": focal["strand"]}
        regions.append(focal_region)
        if candidate["compact_A_B_P-like-1_P-like-2"] != "yes":
            continue
        side = int(candidate["cluster_relative_side"])
        role_specs = (("IdrB_related", side, "IdrB-related"),
                      ("P_like", 2 * side, "P-like-1"), ("P_like", 3 * side, "P-like-2"))
        cluster_regions = [focal_region]
        for family, relative_orf, role in role_specs:
            match = next(row for row in grouped[candidate_id]
                         if row["neighbor_family"] == family
                         and row["neighbor_annotation_accepted"] == "yes"
                         and int(row["relative_ORF"]) == relative_orf)
            region = {**base, "gene_role": role, "gene_id": match["neighbor_id"],
                      "start": int(match["start"]), "end": int(match["end"]), "strand": match["strand"]}
            regions.append(region)
            cluster_regions.append(region)
        clusters[candidate_id] = sorted(cluster_regions, key=lambda item: item["start"])
    return regions, clusters


def write_bed(path, regions):
    with Path(path).open("w") as handle:
        for row in sorted(regions, key=lambda x: (x["contig_id"], x["start"], x["end"])):
            handle.write(f'{row["contig_id"]}\t{row["start"] - 1}\t{row["end"]}\n')


def run_depth(samtools, bam, regions, profile_args):
    depths = defaultdict(dict)
    by_contig = defaultdict(list)
    for row in regions:
        by_contig[row["contig_id"]].append(row)
    for contig, contig_regions in by_contig.items():
        start = min(row["start"] for row in contig_regions)
        end = max(row["end"] for row in contig_regions)
        command = [samtools, "depth", "-a", "-r", f"{contig}:{start}-{end}",
                   *profile_args, str(bam)]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
        for line in process.stdout:
            ref, position, depth = line.rstrip().split("\t")[:3]
            depths[ref][int(position)] = int(depth)
        if process.wait() != 0:
            raise SystemExit(f"samtools depth failed: {' '.join(command)}")
    return depths


def longest_zero(depths):
    longest = current = 0
    for depth in depths:
        if depth == 0:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def metrics(region, depths, profile):
    values = [depths.get(region["contig_id"], {}).get(position, 0)
              for position in range(region["start"], region["end"] + 1)]
    length = len(values)
    third = max(1, length // 3)
    segments = (values[:third], values[third:2 * third], values[2 * third:])
    mean = statistics.fmean(values) if values else 0
    sd = statistics.pstdev(values) if len(values) > 1 else 0
    result = {**region, "profile": profile, "gene_length_bp": length,
              "mean_depth": round(mean, 6), "median_depth": statistics.median(values) if values else 0,
              "max_depth": max(values, default=0), "depth_CV": round(sd / mean, 6) if mean else "",
              "longest_zero_interval_bp": longest_zero(values)}
    for threshold in (1, 2, 5, 10):
        result[f"breadth_{threshold}x"] = round(sum(value >= threshold for value in values) / length, 6) if length else 0
    for label, segment in zip(("five_prime", "middle", "three_prime"), segments):
        result[f"{label}_breadth_1x"] = round(sum(value >= 1 for value in segment) / len(segment), 6) if segment else 0
        result[f"{label}_breadth_5x"] = round(sum(value >= 5 for value in segment) / len(segment), 6) if segment else 0
    return result, values


def reference_span(position, cigar):
    consumed = sum(int(length) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X")
    return position + consumed - 1


def overlap(alignment, region):
    return alignment["start"] <= region["end"] and alignment["end"] >= region["start"]


def junction_evidence(samtools, bam, clusters):
    rows = []
    for candidate_id, genes in clusters.items():
        contig = genes[0]["contig_id"]
        region = f"{contig}:{genes[0]['start']}-{genes[-1]['end']}"
        process = subprocess.run([samtools, "view", "-q", "20", str(bam), region],
                                 check=True, text=True, capture_output=True)
        alignments = defaultdict(list)
        for line in process.stdout.splitlines():
            fields = line.split("\t")
            flag, position, mapq = int(fields[1]), int(fields[3]), int(fields[4])
            if flag & EXCLUDE_PRIMARY:
                continue
            alignment = {"flag": flag, "start": position, "end": reference_span(position, fields[5]),
                         "mapq": mapq, "tlen": abs(int(fields[8]))}
            alignments[fields[0]].append(alignment)
        for left, right in zip(genes, genes[1:]):
            spanning_reads, spanning_pairs, proper_pairs, orientation_pairs = set(), set(), set(), set()
            for read_name, read_alignments in alignments.items():
                if any(aln["start"] <= left["end"] and aln["end"] >= right["start"] for aln in read_alignments):
                    spanning_reads.add(read_name)
                left_hits = [aln for aln in read_alignments if overlap(aln, left)]
                right_hits = [aln for aln in read_alignments if overlap(aln, right)]
                if left_hits and right_hits and len(read_alignments) >= 2:
                    spanning_pairs.add(read_name)
                    if any((a["flag"] & 0x2) and (b["flag"] & 0x2) for a in left_hits for b in right_hits):
                        proper_pairs.add(read_name)
                    if any(bool(a["flag"] & 0x10) != bool(b["flag"] & 0x10) for a in left_hits for b in right_hits):
                        orientation_pairs.add(read_name)
            rows.append({
                "candidate_id": candidate_id, "contig_id": contig,
                "boundary_id": f'{left["gene_role"]}|{right["gene_role"]}',
                "left_gene_id": left["gene_id"], "right_gene_id": right["gene_id"],
                "intergenic_distance_bp": right["start"] - left["end"] - 1,
                "MAPQ20_primary_duplicate_excluded_spanning_reads": len(spanning_reads),
                "MAPQ20_primary_duplicate_excluded_spanning_pairs": len(spanning_pairs),
                "proper_pair_count": len(proper_pairs), "orientation_consistent_pair_count": len(orientation_pairs),
                "analysis_context": "DNA_presence_control_not_transcriptional_evidence",
            })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--neighborhoods", required=True, type=Path)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--sample", default="22_N30_16")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    subprocess.run([args.samtools, "quickcheck", "-v", str(args.bam)], check=True)
    regions, clusters = make_regions(read_tsv(args.classification), read_tsv(args.neighborhoods), args.sample)
    bed = args.outdir / f"{args.sample}_candidate_and_four_gene_regions.bed"
    write_bed(bed, regions)
    all_metrics, mapq20_values = [], {}
    for profile, profile_args in DEPTH_PROFILES.items():
        depths = run_depth(args.samtools, args.bam, regions, profile_args)
        for region in regions:
            row, values = metrics(region, depths, profile)
            all_metrics.append(row)
            if profile == "primary_MAPQ20_deduplicated":
                mapq20_values[(region["candidate_id"], region["gene_role"])] = values
    write_tsv(args.outdir / f"{args.sample}_four_gene_DNA_presence_metrics.tsv", all_metrics)
    with gzip.open(args.outdir / f"{args.sample}_four_gene_per_base_depth_primary_MAPQ20.tsv.gz", "wt") as handle:
        handle.write("candidate_id\tgene_role\tgene_id\tcontig_id\tposition\tdepth\tanalysis_context\n")
        for region in regions:
            values = mapq20_values[(region["candidate_id"], region["gene_role"])]
            for position, depth in zip(range(region["start"], region["end"] + 1), values):
                handle.write(f'{region["candidate_id"]}\t{region["gene_role"]}\t{region["gene_id"]}\t{region["contig_id"]}\t{position}\t{depth}\tDNA_presence_control_not_transcriptional_evidence\n')
    write_tsv(args.outdir / f"{args.sample}_four_gene_junction_DNA_support.tsv",
              junction_evidence(args.samtools, args.bam, clusters))
    write_tsv(args.outdir / f"{args.sample}_four_gene_regions_manifest.tsv", regions)
    report = [
        f"# {args.sample} DNA presence control", "",
        "This analysis uses source metagenomic DNA reads competitively mapped against the whole assembly.",
        "It validates sequence presence and assembly support; it is not transcriptional evidence.", "",
        f"- Candidate IdrA/DMSOR CDS regions: {sum(r['gene_role'] == 'IdrA_candidate' for r in regions)}",
        f"- Complete four-gene neighborhoods evaluated: {len(clusters)}",
        f"- Indexed whole-assembly BAM: `{args.bam}`", "",
    ]
    (args.outdir / f"{args.sample}_DNA_PRESENCE_CONTROL_REPORT.md").write_text("\n".join(report))


if __name__ == "__main__":
    main()
