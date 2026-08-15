#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


FIELDS = [
    "sample", "protein_id", "contig_id", "tree_review_reason", "provisional_tree_class",
    "nearest_anchor_protein", "nearest_anchor_reviewed_class", "nearest_anchor_distance",
    "second_nearest_anchor_class", "second_class_distance_margin", "candidate_anchor_MRCA_support",
    "anchor_classes_below_MRCA", "combined_bitscore", "iriA_bitscore", "aioA_bitscore",
    "combined_HMM_coverage", "protein_length_aa", "prodigal_partial_code",
    "compact_A_B_P-like-1_P-like-2", "IdrB_related_present", "P_like_count", "edge_censored",
    "bin_id", "MAG_quality_class", "GTDB_family", "GTDB_genus", "final_integrated_class",
    "manual_final_class", "manual_reviewer", "manual_evidence_note",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    with args.classification.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t")
                if row["tree_manual_review_required"] == "yes"]
    if len(rows) != 34:
        raise SystemExit(f"Expected 34 pending candidates, found {len(rows)}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "manual_final_class": "pending",
                             "manual_reviewer": "", "manual_evidence_note": ""})
    print(f"pending_candidates={len(rows)}")
    print(f"complete_four_gene_clusters={sum(row['compact_A_B_P-like-1_P-like-2'] == 'yes' for row in rows)}")


if __name__ == "__main__":
    main()
