#!/usr/bin/env python3
"""Audit sample FAA files and build a traceable sample-prefixed protein library."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-faa", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--excluded-file", action="append", default=[])
    return parser.parse_args()


def sample_id(path: Path) -> str:
    suffix = "_AA.faa"
    if not path.name.endswith(suffix):
        raise ValueError(path.name)
    return path.name[: -len(suffix)]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = arguments()
    files = sorted(args.input_dir.glob("*.faa"))
    excluded = set(args.excluded_file)
    if len(files) != 57:
        raise SystemExit(f"Expected 57 FAA files before QC, observed {len(files)}")
    unknown_exclusions = excluded - {path.name for path in files}
    if unknown_exclusions:
        raise SystemExit(f"Excluded files not found: {sorted(unknown_exclusions)}")

    included_samples = [sample_id(path) for path in files if path.name not in excluded]
    if len(included_samples) != len(set(included_samples)):
        raise SystemExit("Duplicate sample IDs remain after exclusions")

    args.output_faa.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    output_digest = hashlib.sha256()
    total_proteins = 0

    with args.output_faa.open("wb") as output:
        for path in files:
            status = "excluded_truncated_duplicate" if path.name in excluded else "included"
            count = 0
            duplicate_headers = 0
            malformed_lines = 0
            seen_ids: set[str] = set()
            saw_header = False
            last_line_ended_newline = True
            with path.open("rb") as handle:
                for raw in handle:
                    last_line_ended_newline = raw.endswith(b"\n")
                    if raw.startswith(b">"):
                        saw_header = True
                        count += 1
                        header = raw[1:].rstrip(b"\r\n")
                        parts = header.split(maxsplit=1)
                        original_id = parts[0].decode("utf-8")
                        if original_id in seen_ids:
                            duplicate_headers += 1
                        else:
                            seen_ids.add(original_id)
                        if status == "included":
                            prefixed = f">{sample_id(path)}|{original_id}".encode()
                            if len(parts) == 2:
                                prefixed += b" " + parts[1]
                            prefixed += b"\n"
                            output.write(prefixed)
                            output_digest.update(prefixed)
                    else:
                        if not saw_header and raw.strip():
                            malformed_lines += 1
                        if status == "included":
                            output.write(raw)
                            output_digest.update(raw)
            if status == "included":
                total_proteins += count
            rows.append({
                "file_name": path.name,
                "sample_id": sample_id(path),
                "status": status,
                "file_size_bytes": path.stat().st_size,
                "protein_count": count,
                "duplicate_header_count_within_file": duplicate_headers,
                "malformed_lines_before_first_header": malformed_lines,
                "ends_with_newline": "yes" if last_line_ended_newline else "no",
                "input_sha256": sha256_file(path),
            })
            print(path.name, status, count, duplicate_headers, flush=True)

    with args.manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.output_faa.parent / "prefixed_contig_all_samples.sha256").write_text(
        f"{output_digest.hexdigest()}  {args.output_faa.name}\n"
    )
    (args.output_faa.parent / "prefixed_contig_counts.tsv").write_text(
        "metric\tvalue\n"
        f"FAA_files_before_QC\t{len(files)}\n"
        f"FAA_files_included\t{len(included_samples)}\n"
        f"unique_samples_included\t{len(set(included_samples))}\n"
        f"total_proteins_included\t{total_proteins}\n"
        f"prefixed_output_size_bytes\t{args.output_faa.stat().st_size}\n"
    )
    if any(row["duplicate_header_count_within_file"] for row in rows if row["status"] == "included"):
        raise SystemExit("Duplicate headers detected in an included FAA file")
    if any(row["ends_with_newline"] == "no" for row in rows if row["status"] == "included"):
        raise SystemExit("An included FAA file is truncated or lacks a final newline")
    print(f"Included {len(included_samples)} samples and {total_proteins} proteins")


if __name__ == "__main__":
    main()
