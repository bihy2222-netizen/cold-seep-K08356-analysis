#!/usr/bin/env python3
"""Create a unique-header whole-assembly reference and reversible crosswalk."""

import argparse
import csv
import re
from collections import Counter


def safe(text): return re.sub(r"[^A-Za-z0-9._:-]+", "_", text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True); parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--crosswalk", required=True)
    args = parser.parse_args(); counts, rows, index = Counter(), [], 0
    out = open(args.out_fasta, "w")
    with open(args.assembly, errors="replace") as source:
        for line in source:
            if line.startswith(">"):
                index += 1; original = line[1:].split()[0]; counts[original] += 1
                renamed = f"22_N30_16|assembly_contig={index:09d}|original={safe(original)}|occurrence={counts[original]}"
                rows.append({"reference_id":renamed,"original_contig_id":original,
                             "original_occurrence":counts[original],"duplicate_original_id":"pending"})
                out.write(f">{renamed}\n")
            else: out.write(line)
    out.close()
    for row in rows: row["duplicate_original_id"] = "yes" if counts[row["original_contig_id"]] > 1 else "no"
    with open(args.crosswalk, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t"); writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__": main()

