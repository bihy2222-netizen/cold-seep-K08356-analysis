#!/usr/bin/env python3
"""Extract depth only for candidate contigs from a whole-assembly BAM."""

import argparse
import csv
import gzip
import subprocess
from collections import defaultdict


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--bam",required=True)
    parser.add_argument("--crosswalk",required=True);parser.add_argument("--strict-genes",required=True)
    parser.add_argument("--out",required=True);args=parser.parse_args()
    cross=defaultdict(list)
    with open(args.crosswalk,newline="") as handle:
        for row in csv.DictReader(handle,delimiter="\t"):cross[row["original_contig_id"]].append(row["reference_id"])
    with open(args.strict_genes,newline="") as handle:
        wanted=sorted({row["contig_id"] for row in csv.DictReader(handle,delimiter="\t")})
    refs=[]
    for contig in wanted:
        if len(cross[contig])!=1:raise SystemExit(f"Whole-assembly reference is not unique for {contig}: {len(cross[contig])}")
        refs.append(cross[contig][0])
    with gzip.open(args.out,"wt") as output:
        for ref in refs:
            process=subprocess.Popen(["samtools","depth","-aa","-r",ref,args.bam],text=True,stdout=subprocess.PIPE)
            for line in process.stdout: output.write(line)
            if process.wait():raise SystemExit(f"samtools depth failed for {ref}")


if __name__=="__main__":main()
