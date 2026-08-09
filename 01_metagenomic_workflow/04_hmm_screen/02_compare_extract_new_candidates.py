#!/usr/bin/env python3
"""Compare full-MAG HMM hits with a known candidate FASTA and extract new hits."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def read_fasta(path: Path, wanted: set[str]) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    chunks: list[str] = []

    def save() -> None:
        if header is None:
            return
        target_id = header.split()[0]
        if target_id in wanted:
            records[target_id] = (header, "".join(chunks))

    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                save()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.strip())
        save()
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--known-fasta", required=True, type=Path)
    parser.add_argument("--combined-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    known = fasta_ids(args.known_fasta)
    by_target: dict[str, dict[str, object]] = {}

    with args.hits.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target = row["prefixed_target_ID"]
            entry = by_target.setdefault(
                target,
                {
                    "MAG_ID": row["MAG_ID"],
                    "protein_ID": row["protein_ID"],
                    "prefixed_target_ID": target,
                    "in_known_candidates": row["protein_ID"] in known,
                    "scores": {},
                },
            )
            scores = entry["scores"]
            assert isinstance(scores, dict)
            scores[row["screen"]] = float(row["full_sequence_bitscore"])

    new_targets = {
        target
        for target, entry in by_target.items()
        if not bool(entry["in_known_candidates"])
    }
    sequences = read_fasta(args.combined_fasta, new_targets)
    missing = sorted(new_targets - sequences.keys())
    if missing:
        raise SystemExit(f"Missing {len(missing)} new targets in combined FASTA: {missing[:5]}")

    screens = [
        "combined_AioA_IdrA_T640",
        "IdrA_specific_T640",
        "AioA_specific_T640",
    ]
    columns = [
        "MAG_ID",
        "protein_ID",
        "prefixed_target_ID",
        "sequence_length_aa",
        "in_known_candidates",
        "evidence_tier",
    ]
    for screen in screens:
        columns.extend([f"{screen}_hit", f"{screen}_bitscore"])

    output_rows: list[dict[str, object]] = []
    for target, entry in sorted(by_target.items()):
        scores = entry["scores"]
        assert isinstance(scores, dict)
        is_known = bool(entry["in_known_candidates"])
        if is_known:
            tier = "existing_candidate"
        elif "combined_AioA_IdrA_T640" in scores:
            tier = "primary_combined_T640_new"
        else:
            tier = "exploratory_IdrA_only_T640_new"
        row: dict[str, object] = {
            "MAG_ID": entry["MAG_ID"],
            "protein_ID": entry["protein_ID"],
            "prefixed_target_ID": target,
            "sequence_length_aa": len(sequences[target][1]) if target in sequences else "",
            "in_known_candidates": is_known,
            "evidence_tier": tier,
        }
        for screen in screens:
            row[f"{screen}_hit"] = screen in scores
            row[f"{screen}_bitscore"] = scores.get(screen, "")
        output_rows.append(row)

    with (args.output_dir / "all_hits_vs_known_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    new_rows = [row for row in output_rows if not bool(row["in_known_candidates"])]
    with (args.output_dir / "new_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(new_rows)

    with (args.output_dir / "new_candidates.faa").open("w") as handle:
        for row in new_rows:
            target = str(row["prefixed_target_ID"])
            header, sequence = sequences[target]
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")

    primary = sum(row["evidence_tier"] == "primary_combined_T640_new" for row in new_rows)
    exploratory = sum(
        row["evidence_tier"] == "exploratory_IdrA_only_T640_new" for row in new_rows
    )
    with (args.output_dir / "comparison_summary.tsv").open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"known_candidate_ids\t{len(known)}\n")
        handle.write(f"union_HMM_targets\t{len(by_target)}\n")
        handle.write(f"new_union_targets\t{len(new_rows)}\n")
        handle.write(f"primary_combined_T640_new\t{primary}\n")
        handle.write(f"exploratory_IdrA_only_T640_new\t{exploratory}\n")


if __name__ == "__main__":
    main()
