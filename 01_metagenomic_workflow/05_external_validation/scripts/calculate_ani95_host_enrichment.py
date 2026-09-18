#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path

from calculate_background_host_enrichment import (
    PARTIAL,
    STRICT,
    enrichment_rows,
    read_tsv,
    write_tsv,
)


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def genome_name_to_key(name):
    stem = Path(name).name
    if stem.endswith(".fa"):
        stem = stem[:-3]
    sample = stem.split(".", 1)[0]
    return f"{sample}|{stem}"


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("-inf")


def representative_sort_key(row):
    return (
        -as_float(row["completeness"]),
        as_float(row["contamination"]),
        -as_float(row["contig_N50"]),
        row["MAG_key"],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--background-taxonomy", required=True, type=Path)
    parser.add_argument("--cdb", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    background = read_tsv(args.background_taxonomy)
    background_by_key = {row["MAG_key"]: row for row in background}
    if len(background_by_key) != 255:
        raise SystemExit(f"Expected 255 background MAGs, found {len(background_by_key)}")

    clusters = defaultdict(list)
    for row in read_csv(args.cdb):
        key = genome_name_to_key(row["genome"])
        if key not in background_by_key:
            raise SystemExit(f"dRep genome not found in background table: {row['genome']} -> {key}")
        clusters[row["secondary_cluster"]].append(background_by_key[key])

    if sum(map(len, clusters.values())) != 255:
        raise SystemExit("ANI95 cluster membership does not contain exactly 255 MAGs")

    membership, representatives = [], []
    for cluster_id, members in sorted(clusters.items()):
        representative = sorted(members, key=representative_sort_key)[0]
        strict_members = [row for row in members if row["strict_IdrA_carrier"] == "yes"]
        partial_members = [row for row in members if row["partial_IdrA_carrier"] == "yes"]
        member_families = sorted({row["GTDB_family"] or "unclassified" for row in members})
        carrier_families = sorted({row["GTDB_family"] or "unclassified"
                                   for row in strict_members + partial_members})

        for row in members:
            membership.append({
                "ANI95_cluster": cluster_id,
                "cluster_size": len(members),
                "sample": row["sample"],
                "MAG_id": row["MAG_id"],
                "MAG_key": row["MAG_key"],
                "completeness": row["completeness"],
                "contamination": row["contamination"],
                "contig_N50": row["contig_N50"],
                "GTDB_family": row["GTDB_family"],
                "strict_IdrA_carrier": row["strict_IdrA_carrier"],
                "partial_IdrA_carrier": row["partial_IdrA_carrier"],
                "selected_representative": "yes" if row is representative else "no",
            })

        representative_row = dict(representative)
        representative_row.update({
            "ANI95_cluster": cluster_id,
            "cluster_size": len(members),
            "cluster_samples": ";".join(sorted({row["sample"] for row in members})),
            "cluster_member_MAG_keys": ";".join(sorted(row["MAG_key"] for row in members)),
            "cluster_strict_carrier": "yes" if strict_members else "no",
            "cluster_partial_carrier": "yes" if partial_members else "no",
            "representative_strict_carrier": representative["strict_IdrA_carrier"],
            "representative_partial_carrier": representative["partial_IdrA_carrier"],
            "cluster_GTDB_families": ";".join(member_families),
            "carrier_GTDB_families": ";".join(carrier_families),
            "family_discordant_within_cluster": "yes" if len(member_families) > 1 else "no",
        })
        # Enrichment is defined at the cluster level, so carrier status follows any member.
        representative_row["strict_IdrA_carrier"] = representative_row["cluster_strict_carrier"]
        representative_row["partial_IdrA_carrier"] = representative_row["cluster_partial_carrier"]
        representatives.append(representative_row)

    write_tsv(args.outdir / "05_ANI95_cluster_membership.tsv", membership)
    write_tsv(args.outdir / "06_ANI95_representative_MAGs.tsv", representatives)

    strict_ids = {row["MAG_key"] for row in representatives if row["strict_IdrA_carrier"] == "yes"}
    extended_ids = {row["MAG_key"] for row in representatives
                    if row["strict_IdrA_carrier"] == "yes" or row["partial_IdrA_carrier"] == "yes"}
    strict = [enrichment_rows(representatives, strict_ids, "strict_ANI95_cluster", unclassified=mode)
              for mode in ("exclude", "include_as_other")]
    extended = [enrichment_rows(representatives, extended_ids,
                                "strict_plus_partial_ANI95_cluster", unclassified=mode)
                for mode in ("exclude", "include_as_other")]
    write_tsv(args.outdir / "07_strict_dereplicated_Fisher_OR_CI.tsv", strict)
    write_tsv(args.outdir / "10_strict_plus_partial_ANI95_sensitivity.tsv", extended)

    summary = [
        {"metric": "input_MQHQ_MAGs", "value": len(background)},
        {"metric": "ANI95_clusters", "value": len(representatives)},
        {"metric": "strict_carrier_clusters", "value": len(strict_ids)},
        {"metric": "strict_plus_partial_carrier_clusters", "value": len(extended_ids)},
        {"metric": "multi_member_clusters", "value": sum(len(v) > 1 for v in clusters.values())},
        {"metric": "family_discordant_clusters", "value": sum(
            row["family_discordant_within_cluster"] == "yes" for row in representatives)},
    ]
    write_tsv(args.outdir / "12_ANI95_analysis_summary.tsv", summary)


if __name__ == "__main__":
    main()
