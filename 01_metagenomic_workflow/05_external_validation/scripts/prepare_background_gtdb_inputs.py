#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 255:
        raise SystemExit(f"Expected 255 MQ/HQ MAGs, found {len(rows)}")

    paths = [row["bin_fasta"] for row in rows]
    sample_mag = [(row["sample"], row["bin_id"]) for row in rows]
    missing = [path for path in paths if not Path(path).is_file()]
    empty = [path for path in paths if Path(path).is_file() and Path(path).stat().st_size == 0]
    if missing or empty or len(set(paths)) != 255 or len(set(sample_mag)) != 255:
        raise SystemExit(json.dumps({"missing": missing, "empty": empty,
                                    "unique_paths": len(set(paths)),
                                    "unique_sample_mag": len(set(sample_mag))}, indent=2))

    crosswalk = []
    for row in rows:
        genome_id = f'{row["sample"]}__{row["bin_id"]}'
        crosswalk.append({**row, "gtdb_user_genome": genome_id})

    batchfile = args.outdir / "all_MQHQ_MAG_GTDB_batchfile.tsv"
    with batchfile.open("w") as handle:
        for row in crosswalk:
            handle.write(f'{row["bin_fasta"]}\t{row["gtdb_user_genome"]}\n')

    fields = list(crosswalk[0])
    crosswalk_path = args.outdir / "all_MQHQ_MAG_GTDB_input_crosswalk.tsv"
    with crosswalk_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(crosswalk)

    with (args.outdir / "drep_genome_paths.txt").open("w") as handle:
        handle.writelines(f'{row["bin_fasta"]}\n' for row in crosswalk)
    with (args.outdir / "drep_genomeInfo.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome", "completeness", "contamination"])
        writer.writeheader()
        for row in crosswalk:
            writer.writerow({"genome": Path(row["bin_fasta"]).name,
                             "completeness": row["Completeness"],
                             "contamination": row["Contamination"]})

    audit = {
        "manifest": str(args.manifest),
        "manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        "MAG_count": len(rows),
        "sample_counts": dict(Counter(row["sample"] for row in rows)),
        "quality_counts": dict(Counter(row["quality_class"] for row in rows)),
        "unique_paths": len(set(paths)),
        "unique_sample_MAG_ids": len(set(sample_mag)),
        "missing_files": len(missing),
        "empty_files": len(empty),
        "total_bytes": sum(Path(path).stat().st_size for path in paths),
        "analysis_unit": "quality-filtered original MetaBAT2 bins; not ANI-dereplicated",
    }
    (args.outdir / "00_input_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
