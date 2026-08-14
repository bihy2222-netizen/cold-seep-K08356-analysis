#!/usr/bin/env python3
"""Plot competitive DMSOR assignments and sparse IdrA-associated coverage."""

import argparse
import csv
import math
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle


SAMPLES = [
    ("SRR16610255_QDN-W01B-1949", "QDN-W01B-1949"),
    ("SRR16610253_QDN-W04B-4900", "QDN-W04B-4900"),
    ("SRR19238834_JL_0.1", "JL_0.1"),
]

FAMILIES = [
    "strict IdrA",
    "partial IdrA-associated",
    "AioA",
    "NapA",
    "ArrA",
    "DmsA",
    "unknown DMSOR",
]

COLORS = {
    "strict IdrA": "#4B136D",
    "partial IdrA-associated": "#A875C3",
    "AioA": "#2F8F5B",
    "NapA": "#E38B2C",
    "ArrA": "#C83E4D",
    "DmsA": "#3478B8",
    "unknown DMSOR": "#8A8A8A",
}

JL_TARGETS = {
    # Sequence-coordinate intervals are from the combined DMSOR HMM alignment
    # in union58_three_model_score_coverage_matrix.tsv (amino acids converted
    # to 1-based nucleotide coordinates).
    "IdrA_PART003": (40, 2718),
    "IdrA_PART021": (28, 2706),
}


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def family_for(reference_id):
    if reference_id.startswith("IdrA_HC"):
        return "strict IdrA"
    if reference_id.startswith("IdrA_PART"):
        return "partial IdrA-associated"
    if "AioA" in reference_id:
        return "AioA"
    if "NapA" in reference_id:
        return "NapA"
    if "ArrA" in reference_id:
        return "ArrA"
    if "DmsA" in reference_id:
        return "DmsA"
    return "unknown DMSOR"


def load_mapping(root):
    mapping_root = root / "maggie_mapping_pilot"
    rows_by_sample = {}
    family_counts = {}
    for sample_id, _ in SAMPLES:
        table = mapping_root / sample_id / "tables" / f"{sample_id}.competitive_63.per_cds_summary.tsv"
        rows = read_tsv(table)
        rows_by_sample[sample_id] = rows
        counts = Counter()
        for row in rows:
            counts[family_for(row["reference_id"])] += int(float(row["mapped_reads"]))
        family_counts[sample_id] = counts
    return rows_by_sample, family_counts


def load_depth(root):
    sample_id = "SRR19238834_JL_0.1"
    path = root / "maggie_mapping_pilot" / sample_id / "tables" / f"{sample_id}.competitive_63.depth.tsv"
    depth = {target: [] for target in JL_TARGETS}
    with path.open() as handle:
        for line in handle:
            reference, position, value = line.rstrip("\n").split("\t")
            if reference in depth:
                depth[reference].append((int(position), int(value)))
    return depth


def load_diamond(root):
    qdn_path = root / "maggie_diamond_blastx_pilot" / "tables" / "diamond_blastx_sample_summary.tsv"
    jl_path = root / "maggie_diamond_blastx_pilot" / "SRR19238834_JL_0.1" / "tables" / "diamond_blastx_sample_summary.tsv"
    combined = {}
    for path in (qdn_path, jl_path):
        for row in read_tsv(path):
            if row.get("read_end") == "combined":
                combined[row["sample"]] = row
    return combined


def write_source_tables(outdir, rows_by_sample, family_counts, depth, diamond):
    family_path = outdir / "panel_A_competitive_family_read_ends.tsv"
    with family_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "family", "mapped_read_ends", "log10_reads_plus_1"])
        for sample_id, label in SAMPLES:
            for family in FAMILIES:
                count = family_counts[sample_id][family]
                writer.writerow([label, family, count, f"{math.log10(count + 1):.6f}"])

    coverage_path = outdir / "panel_B_JL_partial_IdrA_per_base_depth.tsv"
    summary_lookup = {
        row["reference_id"]: row
        for row in rows_by_sample["SRR19238834_JL_0.1"]
        if row["reference_id"] in JL_TARGETS
    }
    with coverage_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["reference_id", "position_nt", "depth", "mapped_read_ends", "coverage_breadth", "hmm_start_nt", "hmm_end_nt"])
        for reference_id, values in depth.items():
            row = summary_lookup[reference_id]
            hmm_start, hmm_end = JL_TARGETS[reference_id]
            for position, value in values:
                writer.writerow([reference_id, position, value, row["mapped_reads"], row["coverage_breadth"], hmm_start, hmm_end])

    evidence_path = outdir / "panel_C_three_sample_evidence_matrix.tsv"
    evidence = build_evidence(rows_by_sample, diamond)
    with evidence_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "strict_reads", "partial_reads", "diamond_support", "max_IdrA_breadth", "partial_read_pairs", "final_classification"])
        for row in evidence:
            writer.writerow([value.replace("\n", " ") for value in row["texts"]])
    return evidence


def build_evidence(rows_by_sample, diamond):
    final_labels = {
        "SRR16610255_QDN-W01B-1949": "other DMSOR\n(reassigned)",
        "SRR16610253_QDN-W04B-4900": "ambiguous\nDMSOR-family",
        "SRR19238834_JL_0.1": "weak partial\n(no strict support)",
    }
    diamond_labels = {
        "SRR16610255_QDN-W01B-1949": "DmsA/NapA; reassigned",
        "SRR16610253_QDN-W04B-4900": "ambiguous DmsA/NapA",
        "SRR19238834_JL_0.1": "ambiguous partial/NapA",
    }
    evidence = []
    for sample_id, label in SAMPLES:
        rows = rows_by_sample[sample_id]
        strict = sum(int(float(row["mapped_reads"])) for row in rows if family_for(row["reference_id"]) == "strict IdrA")
        partial = sum(int(float(row["mapped_reads"])) for row in rows if family_for(row["reference_id"]) == "partial IdrA-associated")
        max_idra = max(
            [float(row["coverage_breadth"]) for row in rows if family_for(row["reference_id"]) in ("strict IdrA", "partial IdrA-associated")] or [0]
        )
        partial_pairs = partial // 2
        scores = [
            min(1.0, strict / 10),
            min(0.65, partial / 12),
            0.15 if "reassigned" in diamond_labels[sample_id] else 0.35,
            min(0.65, max_idra / 0.12),
            min(0.65, partial_pairs / 6),
            0.15 if "reassigned" in final_labels[sample_id] else 0.35,
        ]
        texts = [
            label,
            str(strict),
            str(partial),
            diamond_labels[sample_id],
            f"{100 * max_idra:.1f}%",
            str(partial_pairs),
            final_labels[sample_id],
        ]
        evidence.append({"sample": label, "scores": scores, "texts": texts})
    return evidence


def plot_panel_a(ax, family_counts):
    x = np.arange(len(SAMPLES))
    width = 0.105
    offsets = (np.arange(len(FAMILIES)) - (len(FAMILIES) - 1) / 2) * width
    for offset, family in zip(offsets, FAMILIES):
        values = [math.log10(family_counts[sample_id][family] + 1) for sample_id, _ in SAMPLES]
        ax.bar(x + offset, values, width=width * 0.9, color=COLORS[family], label=family)
    ax.set_xticks(x, [label for _, label in SAMPLES])
    ax.set_ylabel(r"$\log_{10}$(mapped read ends + 1)")
    ax.set_title("A  Competitive nucleotide-level assignments", loc="left", fontweight="bold")
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(ncol=4, frameon=False, fontsize=8, loc="upper left")
    jl_x = x[2]
    for family, count in (("partial IdrA-associated", 6), ("NapA", 2)):
        family_index = FAMILIES.index(family)
        xpos = jl_x + offsets[family_index]
        ypos = math.log10(count + 1)
        ax.text(xpos, ypos + 0.09, str(count), ha="center", va="bottom", fontsize=8, color=COLORS[family], fontweight="bold")
    ax.text(0.995, 0.02, "Bowtie2; counts are mapped read ends", transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#555555")


def plot_depth_axis(ax, reference_id, values, summary_row):
    positions = np.array([item[0] for item in values])
    depths = np.array([item[1] for item in values])
    length = int(summary_row["reference_length"])
    hmm_start, hmm_end = JL_TARGETS[reference_id]
    max_depth = int(depths.max()) if len(depths) else 0
    ax.axvspan(hmm_start, min(hmm_end, length), color="#E38B2C", alpha=0.13, linewidth=0)
    ax.fill_between(positions, depths, step="mid", color=COLORS["partial IdrA-associated"], alpha=0.8)
    ax.plot(positions, depths, color="#693A83", linewidth=0.7)
    ax.add_patch(Rectangle((1, -0.26), length - 1, 0.14, facecolor="#D9D9D9", edgecolor="#777777", linewidth=0.6, clip_on=False))
    covered = depths > 0
    if covered.any():
        starts = positions[covered]
        ax.scatter(starts, np.full(starts.shape, -0.19), marker="s", s=5, color=COLORS["partial IdrA-associated"], clip_on=False)
    ax.set_xlim(1, length)
    ax.set_ylim(-0.3, max(1.3, max_depth * 1.22))
    ax.set_ylabel("Depth")
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    breadth = 100 * float(summary_row["coverage_breadth"])
    read_ends = int(float(summary_row["mapped_reads"]))
    ax.text(
        0.995,
        0.9,
        f"{reference_id}   {read_ends} read ends | breadth {breadth:.1f}% | max depth {max_depth} | MAPQ >= 30",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.3,
    )


def plot_panel_c(ax, evidence):
    matrix = np.array([row["scores"] for row in evidence])
    cmap = LinearSegmentedColormap.from_list("evidence", ["#ECECEC", "#E5A25B", "#C9A7D8", "#4B136D"])
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    columns = ["Strict\nreads", "Partial\nreads", "DIAMOND\nsupport", "Max IdrA\nbreadth", "Partial\nread pairs", "Final call"]
    ax.set_xticks(np.arange(len(columns)), columns)
    ax.set_yticks(np.arange(len(evidence)), [row["sample"] for row in evidence])
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, length=0)
    text_columns = [1, 2, 3, 4, 5, 6]
    for i, row in enumerate(evidence):
        for j, text_index in enumerate(text_columns):
            value = row["texts"][text_index]
            ax.text(j, i, value, ha="center", va="center", fontsize=7.2, color="#222222", wrap=True)
    ax.set_title("C  Integrated evidence strength", loc="left", fontweight="bold", pad=30)
    for edge in np.arange(-0.5, len(columns), 1):
        ax.axvline(edge, color="white", linewidth=1.2)
    for edge in np.arange(-0.5, len(evidence), 1):
        ax.axhline(edge, color="white", linewidth=1.2)
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(evidence), 1), minor=True)
    ax.tick_params(which="minor", bottom=False, left=False)


def make_figure(root, outdir):
    rows_by_sample, family_counts = load_mapping(root)
    depth = load_depth(root)
    diamond = load_diamond(root)
    evidence = write_source_tables(outdir, rows_by_sample, family_counts, depth, diamond)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlepad": 9})
    fig = plt.figure(figsize=(12.4, 12.0), constrained_layout=False)
    grid = fig.add_gridspec(4, 1, height_ratios=[2.7, 1.55, 1.55, 2.25], hspace=0.72)

    ax_a = fig.add_subplot(grid[0])
    plot_panel_a(ax_a, family_counts)

    jl_rows = {row["reference_id"]: row for row in rows_by_sample["SRR19238834_JL_0.1"]}
    ax_b1 = fig.add_subplot(grid[1])
    ax_b1.set_title("B  JL_0.1 partial IdrA-associated per-base coverage", loc="left", fontweight="bold")
    plot_depth_axis(ax_b1, "IdrA_PART003", depth["IdrA_PART003"], jl_rows["IdrA_PART003"])
    ax_b2 = fig.add_subplot(grid[2])
    plot_depth_axis(ax_b2, "IdrA_PART021", depth["IdrA_PART021"], jl_rows["IdrA_PART021"])
    ax_b2.set_xlabel("CDS nucleotide position")
    legend = [
        Patch(facecolor=COLORS["partial IdrA-associated"], label="Observed per-base depth"),
        Patch(facecolor="#E38B2C", alpha=0.25, label="Combined DMSOR HMM-aligned region"),
        Patch(facecolor="#D9D9D9", edgecolor="#777777", label="CDS"),
    ]
    ax_b2.legend(handles=legend, frameon=False, ncol=3, fontsize=8, loc="upper left", bbox_to_anchor=(0, -0.29))

    ax_c = fig.add_subplot(grid[3])
    plot_panel_c(ax_c, evidence)

    fig.suptitle(
        "Competitive metatranscriptomic mapping reveals sparse partial IdrA-associated signals\n"
        "without broad coverage of strict DIRM-like IdrA genes",
        fontsize=14,
        fontweight="bold",
        y=0.988,
    )
    fig.text(
        0.01,
        0.008,
        "Nucleotide counts are mapped read ends; paired support is reported separately. "
        "DIAMOND classified JL_0.1 as ambiguous: partial IdrA had the highest individual score,\n"
        "but NapA dominated assigned read ends (best-family delta 4.6); sparse strict hits did not provide robust strict IdrA support.",
        fontsize=8,
        color="#444444",
    )
    fig.subplots_adjust(top=0.92, bottom=0.07, left=0.12, right=0.985)

    stem = outdir / "competitive_metatranscriptome_IdrA_evidence_three_panel"
    fig.savefig(stem.with_suffix(".png"), dpi=400)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".svg"))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True, help="Maggie expression validation root")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    make_figure(args.root, args.outdir)


if __name__ == "__main__":
    main()
