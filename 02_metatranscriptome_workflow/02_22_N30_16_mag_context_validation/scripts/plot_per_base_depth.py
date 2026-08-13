#!/usr/bin/env python3
"""Plot one SVG per reference from samtools depth -aa output."""

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text)[:180]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--depth", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--title-prefix", default="DNA presence control (not transcriptional evidence)")
    args = parser.parse_args(); out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if args.depth.endswith(".gz") else open
    data = defaultdict(lambda: ([], []))
    with opener(args.depth, "rt") as handle:
        for line in handle:
            ref, position, depth = line.rstrip().split("\t")[:3]
            data[ref][0].append(int(position)); data[ref][1].append(int(depth))
    for ref, (positions, depths) in data.items():
        fig, axis = plt.subplots(figsize=(10, 3.2))
        axis.plot(positions, depths, color="#176B87", linewidth=0.8)
        axis.fill_between(positions, depths, color="#64CCC5", alpha=0.35)
        axis.set(xlabel="Nucleotide position", ylabel="Depth", title=f"{args.title_prefix}: {ref}")
        axis.spines[["top", "right"]].set_visible(False); fig.tight_layout()
        fig.savefig(out / f"{safe_name(ref)}.svg"); fig.savefig(out / f"{safe_name(ref)}.png", dpi=180)
        plt.close(fig)


if __name__ == "__main__": main()
