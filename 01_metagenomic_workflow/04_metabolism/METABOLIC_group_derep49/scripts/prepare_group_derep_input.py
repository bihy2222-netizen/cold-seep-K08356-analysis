#!/usr/bin/env python3
"""Build a checked MAG-protein input set from the four group-derep K08356 hits."""

import argparse
import csv
import hashlib
import os
from collections import defaultdict
from pathlib import Path


GROUPS = {"IS", "AS", "ES", "NS"}
REQUIRED_COLUMNS = {"group", "full_gene_id", "bin_id", "gene_id"}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_count(path):
    with path.open("rb") as handle:
        return sum(1 for line in handle if line.startswith(b">"))


def group_root(binning_root, group):
    return binning_root / group / f"{group}_50_10_dRep" / "dereplicated_genomes"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--binning-root", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expected-mag-count", type=int)
    args = parser.parse_args()

    with args.hits.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Missing required hit-table columns: {sorted(missing)}")
        hits = list(reader)

    if not hits:
        raise SystemExit("The K08356 hit table is empty")
    invalid_groups = sorted({row["group"] for row in hits} - GROUPS)
    if invalid_groups:
        raise SystemExit(f"Unexpected habitat groups: {invalid_groups}")

    by_mag = defaultdict(list)
    for row in hits:
        by_mag[(row["group"], row["bin_id"])].append(row)

    if args.expected_mag_count is not None and len(by_mag) != args.expected_mag_count:
        raise SystemExit(
            f"Expected {args.expected_mag_count} unique MAGs, observed {len(by_mag)}"
        )

    args.input_dir.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    seen_destinations = {}

    for (group, bin_id), mag_hits in sorted(by_mag.items()):
        source = group_root(args.binning_root, group) / f"{bin_id}.cds.faa"
        if not source.is_file() or source.stat().st_size == 0:
            raise SystemExit(f"Missing or empty group-derep protein file: {source}")
        destination = args.input_dir / source.name
        if destination.name in seen_destinations and seen_destinations[destination.name] != source:
            raise SystemExit(f"Duplicate input filename from different sources: {destination.name}")
        seen_destinations[destination.name] = source

        if destination.is_symlink():
            existing = destination.resolve()
            if existing != source.resolve():
                raise SystemExit(f"Existing symlink points elsewhere: {destination} -> {existing}")
        elif destination.exists():
            raise SystemExit(f"Refusing to replace existing non-symlink input: {destination}")
        else:
            os.symlink(source, destination)

        sequence_count = fasta_count(source)
        if sequence_count == 0:
            raise SystemExit(f"No protein sequences found in {source}")
        manifest_rows.append(
            {
                "group": group,
                "bin_id": bin_id,
                "K08356_hit_count": len(mag_hits),
                "K08356_full_gene_ids": ";".join(sorted(row["full_gene_id"] for row in mag_hits)),
                "K08356_gene_ids": ";".join(sorted(row["gene_id"] for row in mag_hits)),
                "source_faa": str(source),
                "input_symlink": str(destination),
                "protein_sequence_count": sequence_count,
                "file_bytes": source.stat().st_size,
                "sha256": sha256(source),
            }
        )

    fields = [
        "group",
        "bin_id",
        "K08356_hit_count",
        "K08356_full_gene_ids",
        "K08356_gene_ids",
        "source_faa",
        "input_symlink",
        "protein_sequence_count",
        "file_bytes",
        "sha256",
    ]
    with args.manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    group_counts = defaultdict(int)
    for row in manifest_rows:
        group_counts[row["group"]] += 1
    print(f"K08356_hits={len(hits)}")
    print(f"unique_MAGs={len(manifest_rows)}")
    print("group_MAG_counts=" + ",".join(f"{group}:{group_counts[group]}" for group in sorted(GROUPS)))
    print(f"input_dir={args.input_dir}")
    print(f"manifest={args.manifest}")


if __name__ == "__main__":
    main()
