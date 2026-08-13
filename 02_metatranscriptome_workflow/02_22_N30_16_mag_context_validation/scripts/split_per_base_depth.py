#!/usr/bin/env python3
"""Split a combined depth -aa stream into one gzip TSV per reference."""

import argparse
import gzip
import re
from pathlib import Path


def safe(text): return re.sub(r"[^A-Za-z0-9._-]+", "_", text)[:200]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--depth", required=True); parser.add_argument("--outdir", required=True)
    args = parser.parse_args(); out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    current, handle = None, None
    with gzip.open(args.depth, "rt") as source:
        for line in source:
            ref = line.split("\t", 1)[0]
            if ref != current:
                if handle: handle.close()
                current = ref; handle = gzip.open(out/f"{safe(ref)}.per_base_depth.tsv.gz", "wt")
                handle.write("reference_id\tposition\tdepth\n")
            handle.write(line)
    if handle: handle.close()


if __name__ == "__main__": main()

