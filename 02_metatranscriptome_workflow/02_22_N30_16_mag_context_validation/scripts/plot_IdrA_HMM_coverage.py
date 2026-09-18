#!/usr/bin/env python3
"""Plot MAPQ-filtered IdrA CDS depth and mark the HMM alignment interval."""

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def safe(text): return re.sub(r"[^A-Za-z0-9._-]+", "_", text)[:180]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--depth", required=True); parser.add_argument("--review", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args(); out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    review = {row["protein_id"]: row for row in csv.DictReader(open(args.review), delimiter="\t")}
    data = defaultdict(lambda: ([], []))
    with gzip.open(args.depth, "rt") as handle:
        for line in handle:
            ref, position, depth = line.rstrip().split("\t")[:3]
            if "|IdrA|" in ref:
                data[ref][0].append(int(position)); data[ref][1].append(int(depth))
    for ref, (positions, depths) in data.items():
        candidate = ref.split("|", 3)[1]; row = review.get(candidate, {})
        start_aa, end_aa = row.get("HMM_alignment_start_aa", ""), row.get("HMM_alignment_end_aa", "")
        fig, axis = plt.subplots(figsize=(10, 3.3))
        axis.plot(positions, depths, color="#8C1D40", linewidth=0.9)
        axis.fill_between(positions, depths, color="#D95F59", alpha=0.3)
        if str(start_aa).isdigit() and str(end_aa).isdigit():
            axis.axvspan((int(start_aa)-1)*3+1, int(end_aa)*3, color="#FFD166", alpha=0.25,
                         label="Combined-HMM alignment interval")
            axis.legend(frameon=False, loc="upper right")
        axis.set(xlabel="IdrA CDS nucleotide position", ylabel="MAPQ>=20 DNA depth",
                 title=f"22_N30_16 DNA presence control (not transcriptional evidence): {candidate}")
        axis.spines[["top","right"]].set_visible(False); fig.tight_layout()
        for ext in ("svg","pdf","png"):
            fig.savefig(out/f"{safe(candidate)}.IdrA_HMM_coverage.{ext}", dpi=180 if ext == "png" else None)
        plt.close(fig)


if __name__ == "__main__": main()
