#!/usr/bin/env python3
"""Summarize alignments and template names separately for each reference."""

import argparse
import csv
import re
import subprocess
from collections import defaultdict


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True); parser.add_argument("--mode", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    data = defaultdict(lambda: {"alignments":0,"read_names":set(),"proper_names":set(),
                                "mapq20_names":set(),"mapq30_names":set(),"nm":0,"aligned":0})
    process = subprocess.Popen(["samtools", "view", args.bam], stdout=subprocess.PIPE, text=True)
    for line in process.stdout:
        fields = line.rstrip().split("\t"); flag, ref, mapq, cigar = int(fields[1]), fields[2], int(fields[4]), fields[5]
        if flag & 4: continue
        stats = data[ref]; stats["alignments"] += 1; stats["read_names"].add(fields[0])
        if flag & 2: stats["proper_names"].add(fields[0])
        if mapq >= 20: stats["mapq20_names"].add(fields[0])
        if mapq >= 30: stats["mapq30_names"].add(fields[0])
        stats["aligned"] += sum(int(n) for n, op in CIGAR.findall(cigar) if op in "M=X")
        for tag in fields[11:]:
            if tag.startswith("NM:i:"): stats["nm"] += int(tag[5:]); break
    if process.wait() != 0: raise SystemExit("samtools view failed")
    fields = ["mode","reference_id","total_alignments","read_names","proper_pair_template_names",
              "MAPQ_ge_20_read_names","MAPQ_ge_30_read_names","mismatches_NM","aligned_bases","mismatch_rate"]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for ref, stats in sorted(data.items()):
            writer.writerow({"mode":args.mode,"reference_id":ref,"total_alignments":stats["alignments"],
                "read_names":len(stats["read_names"]),"proper_pair_template_names":len(stats["proper_names"]),
                "MAPQ_ge_20_read_names":len(stats["mapq20_names"]),"MAPQ_ge_30_read_names":len(stats["mapq30_names"]),
                "mismatches_NM":stats["nm"],"aligned_bases":stats["aligned"],
                "mismatch_rate":f"{stats['nm']/stats['aligned']:.8f}" if stats["aligned"] else "NA"})


if __name__ == "__main__": main()

