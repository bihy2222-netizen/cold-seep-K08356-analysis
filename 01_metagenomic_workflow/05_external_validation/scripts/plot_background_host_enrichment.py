#!/usr/bin/env python3
import argparse
import csv
import math
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_figure(fig, outdir, stem):
    for suffix in ("svg", "pdf", "png"):
        fig.savefig(outdir / f"{stem}.{suffix}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def excluded(rows):
    return [row for row in rows if row["unclassified_handling"] == "exclude"]


def as_float(row, key):
    value = row[key]
    if value.lower() in {"inf", "infinity"}:
        return math.inf
    return float(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", required=True, type=Path)
    parser.add_argument("--ani-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    background = read_tsv(args.bin_dir / "01_all_MQHQ_MAG_GTDB_taxonomy.tsv")
    strict = {row["MAG_key"] for row in background if row["strict_IdrA_carrier"] == "yes"}
    rhodo = {row["MAG_key"] for row in background if row["GTDB_family"] == "f__Rhodobacteraceae"}

    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    groups = [strict, {row["MAG_key"] for row in background} - strict]
    proportions = [100 * len(group & rhodo) / len(group) if group else 0 for group in groups]
    bars = ax.bar(["Strict carriers", "Non-carriers"], proportions,
                  color=["#5B2A86", "#B8B8B8"], width=0.62)
    for bar, group in zip(bars, groups):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.5,
                f"{len(group & rhodo)}/{len(group)}", ha="center", fontsize=9)
    ax.set_ylabel("Rhodobacteraceae (%)")
    ax.set_ylim(0, max(105, max(proportions) + 12))
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, args.outdir, "A_Rhodobacteraceae_proportion")

    family_counts = Counter((row["GTDB_family"] or "unclassified").replace("f__", "")
                            for row in background if row["MAG_key"] in strict)
    labels, counts = zip(*family_counts.most_common()) if family_counts else ([], [])
    fig, ax = plt.subplots(figsize=(6.2, max(3.2, 0.42 * len(labels) + 1.3)))
    y = range(len(labels))
    ax.barh(list(y), counts, color="#5B2A86")
    ax.set_yticks(list(y), labels)
    ax.invert_yaxis()
    ax.set_xlabel("Strict carrier MAGs")
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, args.outdir, "B_strict_carrier_host_families")

    forest_rows = []
    forest_rows += excluded(read_tsv(args.bin_dir / "04_strict_Fisher_OR_CI.tsv"))
    forest_rows += excluded(read_tsv(args.ani_dir / "07_strict_dereplicated_Fisher_OR_CI.tsv"))
    forest_rows += excluded(read_tsv(args.bin_dir / "08_per_sample_enrichment.tsv"))
    forest_rows += excluded(read_tsv(args.bin_dir / "09_leave_one_sample_out_enrichment.tsv"))
    labels = [row["analysis_id"].replace("strict_", "").replace("_", " ") for row in forest_rows]
    odds = [as_float(row, "odds_ratio") for row in forest_rows]
    lows = [as_float(row, "exact_CI95_low") for row in forest_rows]
    highs = [as_float(row, "exact_CI95_high") for row in forest_rows]
    finite_values = [v for v in odds + highs if math.isfinite(v) and v > 0]
    cap = max(10, max(finite_values, default=10) * 1.35)
    plotted_odds = [min(v, cap) if math.isfinite(v) else cap for v in odds]
    plotted_highs = [min(v, cap) if math.isfinite(v) else cap for v in highs]
    fig, ax = plt.subplots(figsize=(8.0, max(4.2, 0.47 * len(labels) + 1.5)))
    y = list(range(len(labels)))
    for i, (point, low, high, original_high) in enumerate(zip(plotted_odds, lows, plotted_highs, highs)):
        low = max(low, 1e-3)
        ax.plot([low, high], [i, i], color="#333333", lw=1.4)
        ax.scatter(point, i, color="#5B2A86", s=34, zorder=3)
        if not math.isfinite(original_high):
            ax.annotate("", xy=(cap, i), xytext=(cap / 1.25, i),
                        arrowprops=dict(arrowstyle="->", color="#333333"))
    ax.axvline(1, color="#888888", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(max(1e-3, min(lows, default=0.1) / 1.4), cap * 1.1)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Odds ratio (exact 95% CI)")
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, args.outdir, "C_strict_OR_forest")

    sensitivity = excluded(read_tsv(args.bin_dir / "04_strict_Fisher_OR_CI.tsv"))
    sensitivity += excluded(read_tsv(args.bin_dir / "10_strict_plus_partial_sensitivity.tsv"))
    sensitivity += excluded(read_tsv(args.ani_dir / "07_strict_dereplicated_Fisher_OR_CI.tsv"))
    sensitivity += excluded(read_tsv(args.ani_dir / "10_strict_plus_partial_ANI95_sensitivity.tsv"))
    labels = [row["analysis_id"].replace("_", " ") for row in sensitivity]
    values = [100 * as_float(row, "Rhodobacteraceae_carrier_prevalence") for row in sensitivity]
    outside = [100 * as_float(row, "outside_Rhodobacteraceae_carrier_prevalence") for row in sensitivity]
    x = list(range(len(labels)))
    fig, ax = plt.subplots(figsize=(8.2, 4.5))
    width = 0.36
    ax.bar([v - width / 2 for v in x], values, width, label="Rhodobacteraceae", color="#5B2A86")
    ax.bar([v + width / 2 for v in x], outside, width, label="Other families", color="#A7A9AC")
    ax.set_xticks(x, labels, rotation=24, ha="right")
    ax.set_ylabel("Carrier prevalence (%)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, args.outdir, "D_strict_partial_sensitivity")


if __name__ == "__main__":
    main()
