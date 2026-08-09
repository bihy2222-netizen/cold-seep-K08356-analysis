#!/usr/bin/env python3
"""Parse EMBOSS needleall output and compare it with the Biopython global run."""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from pathlib import Path


FIELD_RE = re.compile(r"^#\s+([^:]+):\s*(.*)$")
COUNT_RE = re.compile(r"(\d+)\s*/\s*(\d+)\s*\(\s*([0-9.]+)%\)")


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
            row["Sequence_ID"]: row["clade_display"]
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
        row["comparison_category"] = " vs ".join(sorted((str(row["group1"]), str(row["group2"]))))
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
    for category, values in sorted(grouped.items()):
        summary.append({
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
