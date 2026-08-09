#!/usr/bin/env python3
"""Prepare audited MAG/contig candidate FASTAs and fixed small-tree references."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from collections import defaultdict
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag-table", type=Path, required=True)
    parser.add_argument("--mag-faa", type=Path, required=True)
    parser.add_argument("--contig-table", type=Path, required=True)
    parser.add_argument("--contig-faa", type=Path, required=True)
    parser.add_argument("--reference-faa", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def read_fasta(path: Path):
    records = {}
    header = None
    sequence = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if header in records:
                        raise ValueError(f"Duplicate FASTA ID: {header}")
                    records[header] = "".join(sequence).upper()
                header = line[1:].split()[0]
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before first header in {path}")
                sequence.append(line)
    if header is not None:
        if header in records:
            raise ValueError(f"Duplicate FASTA ID: {header}")
        records[header] = "".join(sequence).upper()
    return records


def write_fasta(path: Path, records):
    with path.open("w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def sha256_file(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(row, key):
    value = row.get(key, "")
    if value in ("", "NA", "nan"):
        return float("nan")
    return float(value)


def mag_metadata(table: Path):
    rows = []
    with table.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            scores = [number(row, f"{model}_bitscore") for model in ("IdrA", "AioA", "combined")]
            target = row["prefixed_target_ID"]
            short = f"MAG|{row['MAG_ID']}|{row['protein_ID']}"
            rows.append({
                "source_level": "MAG", "original_target_id": target,
                "short_header": short, "MAG_ID": row["MAG_ID"],
                "sample_id": "", "contig_id": row["protein_ID"].rsplit("_", 1)[0],
                "ORF_ID": row["protein_ID"], "priority": max(scores) >= 640,
                "priority_reason": "at_least_one_exact_model_score_ge640" if max(scores) >= 640 else "T100_archive_only",
                "max_exact_model_score": max(scores),
            })
    return rows


def contig_metadata(table: Path):
    rows = []
    with table.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            scores = [number(row, f"{model}_full_score") for model in ("IdrA", "AioA", "combined")]
            target = row["target_id"]
            short = f"CONTIG|{row['sample_id']}|{row['contig_id']}|{row['orf_id']}"
            rows.append({
                "source_level": "CONTIG", "original_target_id": target,
                "short_header": short, "MAG_ID": "", "sample_id": row["sample_id"],
                "contig_id": row["contig_id"], "ORF_ID": row["orf_id"],
                "priority": max(scores) >= 640,
                "priority_reason": "at_least_one_exact_model_score_ge640" if max(scores) >= 640 else "T100_archive_only",
                "max_exact_model_score": max(scores),
            })
    return rows


def attach_sequences(metadata, fasta):
    missing = []
    output = []
    for row in metadata:
        sequence = fasta.get(row["original_target_id"])
        if sequence is None:
            missing.append(row["original_target_id"])
            continue
        enriched = dict(row)
        enriched["sequence"] = sequence
        enriched["sequence_length"] = len(sequence)
        enriched["sequence_sha256"] = hashlib.sha256(sequence.encode()).hexdigest()
        output.append(enriched)
    return output, missing


def exact_dereplicate(rows):
    clusters = defaultdict(list)
    for row in rows:
        clusters[row["sequence_sha256"]].append(row)
    representatives = []
    membership = []
    for digest, members in sorted(clusters.items()):
        members.sort(key=lambda row: row["short_header"])
        representative = members[0]
        representatives.append((representative["short_header"], representative["sequence"]))
        for member in members:
            membership.append({
                "representative_header": representative["short_header"],
                "member_short_header": member["short_header"],
                "source_level": member["source_level"],
                "MAG_ID": member["MAG_ID"], "sample_id": member["sample_id"],
                "contig_id": member["contig_id"], "ORF_ID": member["ORF_ID"],
                "original_target_id": member["original_target_id"],
                "sequence_sha256": digest,
                "is_representative": "yes" if member is representative else "no",
            })
    representatives.sort()
    return representatives, membership


def write_tsv(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = arguments()
    out = args.outdir.resolve()
    for directory in ("00_audit", "01_NCBI_FASTA", "02_all_T100", "03_derep_tree_input", "04_references"):
        (out / directory).mkdir(parents=True, exist_ok=True)

    mag_fasta = read_fasta(args.mag_faa)
    contig_fasta = read_fasta(args.contig_faa)
    mag_rows, mag_missing = attach_sequences(mag_metadata(args.mag_table), mag_fasta)
    contig_rows, contig_missing = attach_sequences(contig_metadata(args.contig_table), contig_fasta)
    if mag_missing or contig_missing:
        raise SystemExit(f"Missing requested IDs: MAG={len(mag_missing)}, contig={len(contig_missing)}")
    if len(mag_rows) != 720 or len(contig_rows) != 15684:
        raise SystemExit(f"Unexpected candidate counts: MAG={len(mag_rows)}, contig={len(contig_rows)}")

    mag_priority = [row for row in mag_rows if row["priority"]]
    contig_priority = [row for row in contig_rows if row["priority"]]
    mag_priority_reps, mag_priority_members = exact_dereplicate(mag_priority)
    contig_priority_reps, contig_priority_members = exact_dereplicate(contig_priority)

    write_fasta(out / "02_all_T100/MAG_all_T100_candidates.faa", [(row["short_header"], row["sequence"]) for row in mag_rows])
    write_fasta(out / "02_all_T100/contig_all_T100_candidates.faa", [(row["short_header"], row["sequence"]) for row in contig_rows])
    write_fasta(out / "01_NCBI_FASTA/MAG_priority_candidates_all_members.faa", [(row["short_header"], row["sequence"]) for row in mag_priority])
    write_fasta(out / "01_NCBI_FASTA/contig_priority_candidates_all_members.faa", [(row["short_header"], row["sequence"]) for row in contig_priority])
    write_fasta(out / "01_NCBI_FASTA/MAG_priority_candidates_for_NCBI.faa", mag_priority_reps)
    write_fasta(out / "01_NCBI_FASTA/contig_priority_candidates_for_NCBI.faa", contig_priority_reps)

    all_rows = mag_rows + contig_rows
    candidate_reps, candidate_members = exact_dereplicate(all_rows)
    write_fasta(out / "03_derep_tree_input/MAG_contig_candidates_derep.faa", candidate_reps)
    write_tsv(out / "03_derep_tree_input/derep_cluster_members.tsv", candidate_members)
    write_tsv(out / "01_NCBI_FASTA/MAG_priority_header_mapping.tsv", mag_priority_members)
    write_tsv(out / "01_NCBI_FASTA/contig_priority_header_mapping.tsv", contig_priority_members)

    refs = read_fasta(args.reference_faa)
    if len(refs) != 136:
        raise SystemExit(f"Expected 136 fixed small-tree references, observed {len(refs)}")
    ref_copy = out / "04_references/smalltree_reference_sequences.faa"
    shutil.copyfile(args.reference_faa, ref_copy)
    write_fasta(
        out / "03_derep_tree_input/MAG_contig_candidates_derep_with_smalltree_refs.faa",
        candidate_reps + sorted(refs.items()),
    )
    ref_manifest = [
        {"reference_id": header, "length": len(sequence),
         "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()}
        for header, sequence in sorted(refs.items())
    ]
    write_tsv(out / "04_references/smalltree_reference_manifest.tsv", ref_manifest)

    mapping_fields = [
        "source_level", "short_header", "original_target_id", "MAG_ID", "sample_id",
        "contig_id", "ORF_ID", "priority", "priority_reason", "max_exact_model_score",
        "sequence_length", "sequence_sha256",
    ]
    write_tsv(out / "00_audit/header_mapping.tsv", all_rows, mapping_fields)

    output_files = sorted(path for path in out.rglob("*.faa"))
    checksums = [{
        "file": str(path.relative_to(out)), "sequence_count": sum(1 for line in path.open() if line.startswith(">")),
        "sha256": sha256_file(path), "size_bytes": path.stat().st_size,
    } for path in output_files]
    write_tsv(out / "00_audit/output_FASTA_QC_and_SHA256.tsv", checksums)
    summary = [
        {"dataset": "MAG_all_T100", "protein_count": len(mag_rows), "unique_units": len({row['MAG_ID'] for row in mag_rows})},
        {"dataset": "MAG_priority_raw", "protein_count": len(mag_priority), "unique_units": len({row['MAG_ID'] for row in mag_priority})},
        {"dataset": "MAG_priority_NCBI_100pct_derep", "protein_count": len(mag_priority_reps), "unique_units": len({row['MAG_ID'] for row in mag_priority})},
        {"dataset": "contig_all_T100", "protein_count": len(contig_rows), "unique_units": len({row['sample_id'] for row in contig_rows})},
        {"dataset": "contig_priority_raw", "protein_count": len(contig_priority), "unique_units": len({row['sample_id'] for row in contig_priority})},
        {"dataset": "contig_priority_NCBI_100pct_derep", "protein_count": len(contig_priority_reps), "unique_units": len({row['sample_id'] for row in contig_priority})},
        {"dataset": "MAG_contig_all_T100_100pct_derep", "protein_count": len(candidate_reps), "unique_units": "NA"},
        {"dataset": "fixed_smalltree_references", "protein_count": len(refs), "unique_units": "NA"},
    ]
    write_tsv(out / "00_audit/candidate_count_summary.tsv", summary)
    (out / "README_PAUSE_FOR_NCBI.md").write_text(
        "# Candidate preparation pause point\n\n"
        "The MAG and contig T100 candidates were extracted independently. NCBI priority files "
        "use the conservative rule: at least one exact IdrA, AioA, or combined HMM full score "
        ">=640, followed by 100% amino-acid sequence dereplication. All source members remain "
        "in the mapping tables. Broad T100 archives are not IdrA assignments.\n\n"
        "The 136-reference FASTA is copied exactly from the previous small-tree analysis. No "
        "new alignment, tree, functional naming, or gene-neighborhood analysis has been run. "
        "Workflow is paused pending the two NCBI Protein BLAST HitTable files.\n"
    )
    print(summary)


if __name__ == "__main__":
    main()
