#!/usr/bin/env python3
"""Report alignment, template, MAPQ, proper-pair, and mismatch metrics."""

import argparse
import csv
import re
import subprocess


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    alignments = mapq20 = mapq30 = proper_alignments = mismatches = aligned_bases = 0
    primary_names, mapq20_names, proper_names = set(), set(), set()
    process = subprocess.Popen(["samtools", "view", args.bam], stdout=subprocess.PIPE, text=True)
    for line in process.stdout:
        fields = line.rstrip().split("\t"); flag, mapq, cigar = int(fields[1]), int(fields[4]), fields[5]
        if flag & 4: continue
        alignments += 1
        if not flag & (256 | 2048): primary_names.add(fields[0])
        if mapq >= 20: mapq20 += 1; mapq20_names.add(fields[0])
        if mapq >= 30: mapq30 += 1
        if flag & 2: proper_alignments += 1; proper_names.add(fields[0])
        aligned_bases += sum(int(n) for n, op in CIGAR.findall(cigar) if op in "M=X")
        for tag in fields[11:]:
            if tag.startswith("NM:i:"): mismatches += int(tag[5:]); break
    if process.wait() != 0: raise SystemExit("samtools view failed")
    row = {"mode": args.mode, "total_alignments": alignments, "primary_read_names": len(primary_names),
           "MAPQ_ge_20_alignments": mapq20, "MAPQ_ge_20_read_names": len(mapq20_names),
           "MAPQ_ge_30_alignments": mapq30, "proper_pair_alignments": proper_alignments,
           "proper_pair_template_names": len(proper_names), "mismatches_NM": mismatches,
           "aligned_bases": aligned_bases,
           "mismatch_rate": f"{mismatches/aligned_bases:.8f}" if aligned_bases else "NA"}
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row, delimiter="\t"); writer.writeheader(); writer.writerow(row)


if __name__ == "__main__": main()
