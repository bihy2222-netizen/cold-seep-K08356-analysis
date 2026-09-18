#!/usr/bin/env python3
"""Stream samtools depth -aa summaries while preserving zero positions."""

import argparse
import csv
import gzip
import math
import statistics
import sys


FIELDS = ["mode","reference_id","length","covered_bases_ge_1x","covered_bases_ge_2x",
          "covered_bases_ge_5x","covered_bases_ge_10x","breadth_ge_1x","breadth_ge_2x",
          "breadth_ge_5x","breadth_ge_10x","mean_depth","median_depth","maximum_depth",
          "depth_CV","longest_zero_coverage_interval","five_prime_breadth_1x",
          "middle_breadth_1x","three_prime_breadth_1x"]


def open_text(path):
    if str(path) == "-": return sys.stdin
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def summarize(mode, ref, depths):
    length=len(depths);mean=statistics.fmean(depths) if depths else 0;zero=longest=0
    for value in depths:zero=zero+1 if value==0 else 0;longest=max(longest,zero)
    counts={n:sum(value>=n for value in depths) for n in (1,2,5,10)}
    variance=statistics.pvariance(depths) if depths else 0
    a,b=math.ceil(length/3),math.ceil(2*length/3);thirds=(depths[:a],depths[a:b],depths[b:])
    breadths=[sum(value>=1 for value in part)/len(part) if part else 0 for part in thirds]
    return {"mode":mode,"reference_id":ref,"length":length,
        **{f"covered_bases_ge_{n}x":counts[n] for n in (1,2,5,10)},
        **{f"breadth_ge_{n}x":f"{counts[n]/length:.6f}" if length else "0" for n in (1,2,5,10)},
        "mean_depth":f"{mean:.6f}","median_depth":statistics.median(depths) if depths else 0,
        "maximum_depth":max(depths) if depths else 0,"depth_CV":f"{math.sqrt(variance)/mean:.6f}" if mean else "NA",
        "longest_zero_coverage_interval":longest,"five_prime_breadth_1x":f"{breadths[0]:.6f}",
        "middle_breadth_1x":f"{breadths[1]:.6f}","three_prime_breadth_1x":f"{breadths[2]:.6f}"}


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--depth",required=True)
    parser.add_argument("--out",required=True);parser.add_argument("--mode",required=True);args=parser.parse_args()
    current,depths=None,[]
    with open(args.out,"w",newline="") as output:
        writer=csv.DictWriter(output,fieldnames=FIELDS,delimiter="\t");writer.writeheader()
        with open_text(args.depth) as source:
            for line in source:
                ref,_position,depth=line.rstrip().split("\t")[:3]
                if current is not None and ref!=current:
                    writer.writerow(summarize(args.mode,current,depths));depths=[]
                current=ref;depths.append(int(depth))
        if current is not None:writer.writerow(summarize(args.mode,current,depths))


if __name__=="__main__":main()

