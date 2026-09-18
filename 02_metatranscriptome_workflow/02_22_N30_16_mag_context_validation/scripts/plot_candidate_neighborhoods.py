#!/usr/bin/env python3
"""Draw compact SVG/PDF/PNG candidate-neighborhood arrows from reviewed rows."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow


COLORS = {"IdrA":"#B2182B", "IdrB_related":"#2166AC", "P_like":"#1B7837", "other":"#BDBDBD"}


def safe(text): return re.sub(r"[^A-Za-z0-9._-]+", "_", text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--neighborhood", required=True); parser.add_argument("--outdir", required=True)
    args = parser.parse_args(); out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    groups = defaultdict(list)
    with open(args.neighborhood, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"): groups[row["candidate_id"]].append(row)
    for candidate, rows in sorted(groups.items()):
        rows.sort(key=lambda row: int(row["start"])); left = min(int(row["start"]) for row in rows)
        right = max(int(row["end"]) for row in rows); span = max(right-left+1, 1)
        fig, axis = plt.subplots(figsize=(12, 2.8))
        for row in rows:
            family = "IdrA" if row["neighbor_id"] == candidate else row.get("neighbor_family", "other")
            family = family if row.get("neighbor_annotation_accepted") == "yes" or family == "IdrA" else "other"
            start, end = int(row["start"])-left, int(row["end"])-left; width = end-start+1
            if row["strand"] == "+": x, dx = start, width
            else: x, dx = end, -width
            axis.add_patch(FancyArrow(x, 0, dx, 0, width=0.36, head_width=0.58,
                                      head_length=min(max(span*0.012, 20), abs(dx)*0.45),
                                      length_includes_head=True, color=COLORS.get(family, COLORS["other"])))
            label = "IdrA" if family == "IdrA" else family.replace("_related", "").replace("_like", "")
            axis.text((start+end)/2, 0.43, label, ha="center", va="bottom", fontsize=8, rotation=30)
        axis.set_xlim(-span*0.02, span*1.02); axis.set_ylim(-0.7, 1.1); axis.set_yticks([])
        axis.set_xlabel("Position relative to extracted neighborhood (bp)")
        axis.set_title(f"22_N30_16 candidate neighborhood: {candidate}")
        axis.spines[["left","right","top"]].set_visible(False); fig.tight_layout()
        for ext in ("svg","pdf","png"):
            fig.savefig(out/f"{safe(candidate)}.neighborhood.{ext}", dpi=180 if ext == "png" else None)
        plt.close(fig)


if __name__ == "__main__": main()
