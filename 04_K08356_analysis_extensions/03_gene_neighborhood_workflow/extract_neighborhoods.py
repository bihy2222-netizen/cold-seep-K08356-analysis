#!/usr/bin/env python3
import argparse, csv
from collections import defaultdict
from pathlib import Path

def fasta_ids(path):
    for line in Path(path).open(errors="ignore"):
        if line.startswith(">"):
            yield line[1:].split()[0]

def attrs(text):
    out = {}
    for item in text.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v.strip('"')
    return out

ap = argparse.ArgumentParser()
ap.add_argument("--annotation-root", required=True)
ap.add_argument("--targets", required=True)
ap.add_argument("--out-long", required=True)
ap.add_argument("--out-summary", required=True)
args = ap.parse_args()

root = Path(args.annotation_root)
targets = set(fasta_ids(args.targets))
by_id, by_contig = {}, defaultdict(list)

for group in ("IS", "AS", "ES", "NS"):
    for gff in sorted((root / group / "01_prodigal" / "GFF").glob("*.gff")):
        genome = gff.stem
        for line in gff.open(errors="ignore"):
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip().split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            a = attrs(p[8])
            gid = a.get("ID")
            if not gid:
                continue
            row = dict(group=group, genome_id=genome, contig_id=p[0],
                       gene_id=gid, start=int(p[3]), end=int(p[4]),
                       strand=p[6], partial=a.get("partial", ""))
            by_id[gid] = row
            by_contig[(group, genome, p[0])].append(row)

for genes in by_contig.values():
    genes.sort(key=lambda x: (x["start"], x["end"], x["gene_id"]))

fields = ["target_id","group","genome_id","contig_id","relative_position","role",
          "neighbor_raw_id","start","end","strand","partial","left_orfs_available",
          "right_orfs_available","target_to_left_cds_bp","target_to_right_cds_bp",
          "edge_censoring"]
summary_fields = ["target_id","group","genome_id","contig_id","target_found",
                  "left_orfs_available","right_orfs_available",
                  "target_to_left_cds_bp","target_to_right_cds_bp","edge_censoring"]
long_rows, summaries = [], []

for tid in sorted(targets):
    t = by_id.get(tid)
    if not t:
        summaries.append(dict(target_id=tid, target_found="no", edge_censoring="unknown"))
        continue
    genes = by_contig[(t["group"], t["genome_id"], t["contig_id"])]
    idx = next(i for i, g in enumerate(genes) if g["gene_id"] == tid)
    left_n, right_n = idx, len(genes) - idx - 1
    contig_min = min(g["start"] for g in genes)
    contig_max = max(g["end"] for g in genes)
    left_bp = t["start"] - contig_min
    right_bp = contig_max - t["end"]
    left_censored = left_n < 10 or left_bp < 10000
    right_censored = right_n < 10 or right_bp < 10000
    if left_censored and right_censored: censor = "both_edge_censored"
    elif left_censored: censor = "left_edge_censored"
    elif right_censored: censor = "right_edge_censored"
    else: censor = "complete_observable"
    common = dict(target_id=tid, group=t["group"], genome_id=t["genome_id"],
                  contig_id=t["contig_id"], left_orfs_available=left_n,
                  right_orfs_available=right_n, target_to_left_cds_bp=left_bp,
                  target_to_right_cds_bp=right_bp, edge_censoring=censor)
    summaries.append(dict(common, target_found="yes"))
    for j in range(max(0, idx-10), min(len(genes), idx+11)):
        g = genes[j]
        long_rows.append(dict(common, relative_position=j-idx,
                              role="target" if j == idx else "neighbor",
                              neighbor_raw_id=g["gene_id"], start=g["start"],
                              end=g["end"], strand=g["strand"], partial=g["partial"]))

for path, rows, flds in [(args.out_long,long_rows,fields),
                         (args.out_summary,summaries,summary_fields)]:
    with open(path, "w", newline="") as h:
        w = csv.DictWriter(h, delimiter="\t", fieldnames=flds, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
