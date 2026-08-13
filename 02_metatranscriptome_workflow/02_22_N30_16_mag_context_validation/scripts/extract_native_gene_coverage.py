#!/usr/bin/env python3
"""Extract four-gene native coverage from whole-assembly competitive depth."""

import argparse
import csv
import gzip
import math
import statistics
from collections import defaultdict


def read_tsv(path):
    with open(path, newline="") as handle: return list(csv.DictReader(handle, delimiter="\t"))


def metrics(depths):
    length = len(depths); mean = statistics.fmean(depths) if depths else 0
    counts = {n: sum(value >= n for value in depths) for n in (1,2,5,10)}
    zero = longest = 0
    for value in depths: zero = zero + 1 if value == 0 else 0; longest = max(longest, zero)
    a, b = math.ceil(length/3), math.ceil(2*length/3); thirds=(depths[:a],depths[a:b],depths[b:])
    return {**{f"breadth_{n}x":counts[n]/length if length else 0 for n in (1,2,5,10)},
            "mean_depth":mean,"median_depth":statistics.median(depths) if depths else 0,
            "max_depth":max(depths) if depths else 0,
            "depth_CV":statistics.pstdev(depths)/mean if mean else "NA",
            "longest_zero_interval":longest,
            "five_prime_breadth":sum(x>=1 for x in thirds[0])/len(thirds[0]) if thirds[0] else 0,
            "middle_breadth":sum(x>=1 for x in thirds[1])/len(thirds[1]) if thirds[1] else 0,
            "three_prime_breadth":sum(x>=1 for x in thirds[2])/len(thirds[2]) if thirds[2] else 0}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--depth",required=True)
    parser.add_argument("--crosswalk",required=True); parser.add_argument("--strict-genes",required=True)
    parser.add_argument("--out",required=True); args=parser.parse_args()
    cross=defaultdict(list)
    for row in read_tsv(args.crosswalk): cross[row["original_contig_id"]].append(row["reference_id"])
    genes=read_tsv(args.strict_genes); wanted={row["contig_id"] for row in genes}; depths={}
    for contig in wanted:
        refs=cross.get(contig,[])
        if len(refs)!=1: raise SystemExit(f"Whole-assembly contig mapping is not unique for {contig}: {len(refs)}")
        depths[refs[0]]=[]
    with gzip.open(args.depth,"rt") as handle:
        for line in handle:
            ref,_pos,dep=line.rstrip().split("\t")[:3]
            if ref in depths: depths[ref].append(int(dep))
    rows=[]
    for gene in genes:
        ref=cross[gene["contig_id"]][0]; segment=depths[ref][int(gene["start"])-1:int(gene["end"])]
        rows.append({"evidence_level":"A_whole_assembly_primary","sample":"22_N30_16",
                     "candidate_id":gene["candidate_id"],"role":gene["role"],"gene_id":gene["gene_id"],
                     "contig_id":gene["contig_id"],"start":gene["start"],"end":gene["end"],
                     "strand":gene["strand"],**metrics(segment),
                     "interpretation":"DNA_presence_control_not_transcriptional_evidence"})
    fields=list(rows[0]) if rows else []
    with open(args.out,"w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t");writer.writeheader();writer.writerows(rows)


if __name__=="__main__":main()
