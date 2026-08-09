#!/usr/bin/env python3
"""Parse EMBOSS needleall output and compare it with the Biopython global run."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
import statistics
from pathlib import Path


FIELD_RE = re.compile(r"^#\s+([^:]+):\s*(.*)$")
COUNT_RE = re.compile(r"(\d+)\s*/\s*(\d+)\s*\(\s*([0-9.]+)%\)")
GROUP_NAMES = {
    "Canonical AioA": "canonical AioA-associated",
    "DIRM-synteny IdrA": "synteny-supported strict DIRM-like IdrA",
    "IdrA phylogenetic": "partial IdrA-associated",
    "Uncertain DMSOR": "AioA-like or unresolved DMSOR",
}
GROUP_ORDER = list(GROUP_NAMES.values())
PAIR_ORDER = (
    [(group, group) for group in GROUP_ORDER]
    + list(itertools.combinations(GROUP_ORDER, 2))
)


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--needle", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--biopython", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def clean_id(value: str) -> str:
    value = value.strip()
    return value.split("|", 1)[-1]


def fasta_ids(path: Path) -> list[str]:
    ids = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(clean_id(line[1:].split()[0]))
    return ids


def restore_emboss_id(value: str, ids: list[str]) -> str:
    match = re.fullmatch(r"EMBOSS_(\d+)", value)
    if not match:
        return clean_id(value)
    number = int(match.group(1))
    index = (number - 1) % len(ids)
    return ids[index]


def parse_count(value: str) -> tuple[int, int]:
    match = COUNT_RE.search(value)
    if not match:
        raise ValueError(f"Cannot parse EMBOSS count field: {value}")
    return int(match.group(1)), int(match.group(2))


def parse_needle(path: Path, ids: list[str]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    current: dict[str, str] = {}
    with path.open(errors="replace") as handle:
        for raw in handle:
            match = FIELD_RE.match(raw.rstrip("\n"))
            if not match:
                continue
            key, value = match.groups()
            if key == "1" and current:
                if {"1", "2", "Length", "Identity", "Similarity", "Gaps", "Score"} <= current.keys():
                    records.append(convert_record(current, ids))
                current = {}
            if key in {"1", "2", "Length", "Identity", "Similarity", "Gaps", "Score"}:
                current[key] = value
    if current and {"1", "2", "Length", "Identity", "Similarity", "Gaps", "Score"} <= current.keys():
        records.append(convert_record(current, ids))
    return records


def convert_record(row: dict[str, str], ids: list[str]) -> dict[str, object]:
    identical, identity_denominator = parse_count(row["Identity"])
    similar, similarity_denominator = parse_count(row["Similarity"])
    gaps, gap_denominator = parse_count(row["Gaps"])
    length = int(row["Length"])
    if len({length, identity_denominator, similarity_denominator, gap_denominator}) != 1:
        raise ValueError(f"Inconsistent EMBOSS alignment lengths: {row}")
    return {
        "seq1": restore_emboss_id(row["1"], ids),
        "seq2": restore_emboss_id(row["2"], ids),
        "alignment_length": length,
        "identical_count": identical,
        "similar_count": similar,
        "gap_count": gaps,
        "needle_identity_pct": 100.0 * identical / length,
        "needle_similarity_pct": 100.0 * similar / length,
        "needle_gap_pct": 100.0 * gaps / length,
        "needle_score": float(row["Score"]),
    }


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def main() -> None:
    options = args()
    options.outdir.mkdir(parents=True, exist_ok=True)
    with options.metadata.open() as handle:
        metadata = {
            row["Sequence_ID"]: GROUP_NAMES[row["clade_display"]]
            for row in csv.DictReader(handle, delimiter="\t")
        }

    ordered_ids = fasta_ids(options.fasta)
    if len(ordered_ids) != 49 or len(set(ordered_ids)) != 49:
        raise SystemExit(f"Expected 49 unique FASTA IDs, observed {len(ordered_ids)}")
    directed = parse_needle(options.needle, ordered_ids)
    expected_ids = set(metadata)
    seen: dict[tuple[str, str], dict[str, object]] = {}
    for row in directed:
        if row["seq1"] not in expected_ids or row["seq2"] not in expected_ids:
            raise SystemExit(f"Unknown sequence ID in needleall output: {row['seq1']} / {row['seq2']}")
        if row["seq1"] == row["seq2"]:
            continue
        key = tuple(sorted((str(row["seq1"]), str(row["seq2"]))))
        if key in seen:
            previous = seen[key]
            previous["reverse_needle_identity_pct"] = row["needle_identity_pct"]
            previous["directional_identity_difference_pct_points"] = (
                float(row["needle_identity_pct"]) - float(previous["needle_identity_pct"])
            )
            continue
        row["group1"] = metadata[str(row["seq1"])]
        row["group2"] = metadata[str(row["seq2"])]
        ordered_groups = sorted(
            (str(row["group1"]), str(row["group2"]),), key=GROUP_ORDER.index
        )
        row["comparison_category"] = " vs ".join(ordered_groups)
        row["reverse_needle_identity_pct"] = math.nan
        row["directional_identity_difference_pct_points"] = math.nan
        seen[key] = row

    rows = list(seen.values())
    if len(rows) != 1176:
        raise SystemExit(f"Expected 1176 nonredundant pairs, observed {len(rows)}")

    long_path = options.outdir / "EMBOSS_Needle_pairwise_identity_long.tsv"
    fields = list(rows[0])
    with long_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    matrix = {sid: {other: "" for other in metadata} for sid in metadata}
    for sid in metadata:
        matrix[sid][sid] = 100.0
    for row in rows:
        matrix[str(row["seq1"])][str(row["seq2"])] = row["needle_identity_pct"]
        matrix[str(row["seq2"])][str(row["seq1"])] = row["needle_identity_pct"]
    with (options.outdir / "EMBOSS_Needle_identity_matrix.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        ids = list(metadata)
        writer.writerow(["Sequence_ID", *ids])
        for sid in ids:
            writer.writerow([sid, *(matrix[sid][other] for other in ids)])

    grouped: dict[str, list[float]] = {}
    for row in rows:
        grouped.setdefault(str(row["comparison_category"]), []).append(float(row["needle_identity_pct"]))
    summary = []
    for group1, group2 in PAIR_ORDER:
        category = f"{group1} vs {group2}"
        values = grouped[category]
        summary.append({
            "comparison_type": "within-clade" if group1 == group2 else "between-clade",
            "comparison_category": category,
            "n_pairs": len(values),
            "minimum": min(values),
            "P05": quantile(values, 0.05),
            "Q1": quantile(values, 0.25),
            "median": statistics.median(values),
            "mean": statistics.mean(values),
            "Q3": quantile(values, 0.75),
            "P95": quantile(values, 0.95),
            "maximum": max(values),
            "SD": statistics.stdev(values) if len(values) > 1 else 0.0,
        })
    with (options.outdir / "EMBOSS_Needle_clade_identity_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    report_dir = options.outdir.parent / "09_report"
    report_dir.mkdir(parents=True, exist_ok=True)
    report_table = report_dir / "four_clade_complete_EMBOSS_Needle_identity_10_comparisons.tsv"
    with report_table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    lookup = {row["comparison_category"]: row for row in summary}
    canonical = GROUP_ORDER[0]
    strict = GROUP_ORDER[1]
    strict_within = lookup[f"{strict} vs {strict}"]
    canonical_strict = lookup[f"{canonical} vs {strict}"]
    nonoverlap = float(strict_within["minimum"]) - float(canonical_strict["maximum"])
    with (report_dir / "README_four_clade_identity_final_CN.md").open("w") as handle:
        handle.write("# 49 条核心蛋白四 clade 全局 identity 最终结果\n\n")
        handle.write("49 条核心蛋白已依据系统树和基因邻域分为四个 clade；")
        handle.write("本分析不是重新分组。720 条 T100 宽松 HMM 候选不纳入本表。\n\n")
        handle.write("## 固定名称与规模\n\n")
        for group, count in zip(GROUP_ORDER, (3, 20, 19, 7)):
            handle.write(f"- `{group}`: {count}\n")
        handle.write("\n后两组为分析类别，并非实验验证的酶功能分类。\n\n")
        handle.write("## 完整 10 类比较\n\n")
        handle.write("| 类型 | 比较 | n | minimum | Q1 | median | Q3 | maximum | mean | SD |\n")
        handle.write("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in summary:
            handle.write(
                f"| {row['comparison_type']} | {row['comparison_category']} | {row['n_pairs']} | "
                f"{row['minimum']:.2f}% | {row['Q1']:.2f}% | {row['median']:.2f}% | "
                f"{row['Q3']:.2f}% | {row['maximum']:.2f}% | {row['mean']:.2f}% | {row['SD']:.2f} |\n"
            )
        handle.write("\n## 解释边界\n\n")
        handle.write(
            f"strict IdrA 组内最低 identity（{strict_within['minimum']:.2f}%）与 "
            f"AioA-strict IdrA 组间最高 identity（{canonical_strict['maximum']:.2f}%）之间"
            f"存在 {nonoverlap:.2f} 个百分点的无重叠间隔。该数值描述两个分布边界的距离，"
            "不是分类阈值区间。partial 和 unresolved clade 的完整分布用于评估两个主要"
            "端点之间的过渡性与内部异质性。通用功能阈值仍需外部实验验证参考、系统树、"
            "HMM 竞争分值和基因邻域共同验证。\n"
        )

    with options.biopython.open() as handle:
        biopython = {
            tuple(sorted((row["seq1"], row["seq2"]))): float(row["global_identity_with_gaps"])
            for row in csv.DictReader(handle, delimiter="\t")
        }
    comparison = []
    for key, row in seen.items():
        python_value = biopython[key]
        needle_value = float(row["needle_identity_pct"])
        comparison.append({
            "seq1": key[0], "seq2": key[1],
            "EMBOSS_Needle_identity_pct": needle_value,
            "Biopython_global_identity_pct": python_value,
            "Needle_minus_Biopython_pct_points": needle_value - python_value,
        })
    with (options.outdir / "EMBOSS_vs_Biopython_global_identity.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(comparison[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(comparison)
    max_difference = max(abs(float(row["Needle_minus_Biopython_pct_points"])) for row in comparison)
    max_directional_difference = max(
        abs(float(row["directional_identity_difference_pct_points"])) for row in rows
    )
    print(f"Parsed {len(directed)} directed records and {len(rows)} nonredundant pairs")
    print(f"Maximum absolute EMBOSS-vs-Biopython identity difference: {max_difference:.12g} percentage points")
    print(f"Maximum EMBOSS forward-vs-reverse identity difference: {max_directional_difference:.12g} percentage points")


if __name__ == "__main__":
    main()
