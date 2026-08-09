#!/usr/bin/env python3
"""Add the new group-derep K08356 MAG to the existing 48-MAG evidence tables."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE_48 = ROOT / "source_48"
SOURCE_NEW = ROOT / "source_new_MAG"
OUTPUT = ROOT / "outputs"

TARGET = "R2111_N500_0-10_bin13-k141_4404301_3"
TREE_LABEL = f"NS|{TARGET}"
BIN_ID = "R2111_N500_0-10_bin13"
CONTIG_ID = "R2111_N500_0-10_bin13-k141_4404301"
SISTER = "IS|SY368YW-8-12_bin16-k141_314641_4"
NEAREST_REFERENCE = "ISMEJFig5|Unknown_clade|MBT97654.1|Dehalococcoidia_bacterium"


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return rows, list(reader.fieldnames or [])


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def parse_gff(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] != CONTIG_ID or fields[2] != "CDS":
                continue
            attrs = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = value
            prod_id = attrs.get("ID", "")
            if not prod_id:
                continue
            gene_index = prod_id.split("_")[-1]
            gene_id = f"{CONTIG_ID}_{gene_index}"
            records[gene_id] = {
                "gene_id": gene_id,
                "gene_index": gene_index,
                "start": fields[3],
                "end": fields[4],
                "strand": fields[6],
                "partial": attrs.get("partial", ""),
            }
    return records


def parse_kofam(path: Path) -> dict[str, list[dict[str, str]]]:
    results: dict[str, list[dict[str, str]]] = {}
    with path.open(newline="") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if not fields or fields[0].startswith("#") or len(fields) < 7:
                continue
            row = {
                "reliable": "yes" if fields[0] == "*" else "no",
                "gene_id": fields[1],
                "ko": fields[2],
                "threshold": fields[3],
                "score": fields[4],
                "evalue": fields[5],
                "definition": fields[6].strip('"'),
            }
            results.setdefault(row["gene_id"], []).append(row)
    return results


def best_annotation(gene_id: str, kofam: dict[str, list[dict[str, str]]]) -> tuple[str, str]:
    rows = kofam.get(gene_id, [])
    reliable = [row for row in rows if row["reliable"] == "yes"]
    if reliable:
        row = max(reliable, key=lambda item: float(item["score"]))
        annotation = (
            f"rank: KOfam reliable; {row['definition']} "
            f"({row['ko']}; score={row['score']}; threshold={row['threshold']}; e={row['evalue']})"
        )
        label = "aioB" if row["ko"] == "K08355" else ""
        return annotation, label
    if rows:
        row = max(rows, key=lambda item: float(item["score"]))
        return (
            f"unresolved/partial CDS; weak KOfam match below adaptive threshold: "
            f"{row['definition']} ({row['ko']}; score={row['score']}; threshold={row['threshold']})",
            "",
        )
    return "unresolved/partial CDS; no KOfam match", ""


def ref_family(subject_id: str) -> str:
    if subject_id.startswith(("IdrP_like", "IdrP1", "IdrP2")):
        return "P_like"
    if subject_id.startswith("AioB"):
        return "canonical_AioB"
    if subject_id.startswith("IdrB"):
        return "IdrB_related"
    return "other"


def accepted_hit(row: dict[str, str]) -> bool:
    family = row["ref_family"]
    evalue = float(row["evalue"])
    bitscore = float(row["bitscore"])
    qcov = float(row["qcov_pct"])
    scov = float(row["scov_pct"])
    if family in {"IdrB_related", "canonical_AioB"}:
        return evalue <= 1e-5 and bitscore >= 50 and qcov >= 50 and scov >= 35
    if family == "P_like":
        return evalue <= 1e-10 and bitscore >= 80 and qcov >= 45 and scov >= 35
    return False


def parse_new_diamond_hits(annotation: str) -> list[dict[str, str]]:
    fields = [
        "qseqid", "sseqid", "pident", "alignment_length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
    ]
    output = []
    path = SOURCE_NEW / "joint_reference_hits.raw.tsv"
    with path.open(newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if len(values) < 14:
                continue
            raw = dict(zip(fields, values[:14]))
            qlen = float(raw["qlen"])
            slen = float(raw["slen"])
            length = float(raw["alignment_length"])
            row = {
                **raw,
                "target_id": TARGET,
                "bin_id": BIN_ID,
                "contig_id": CONTIG_ID,
                "relative_position": "-1",
                "role": "upstream",
                "neighbor_raw_id": raw["qseqid"],
                "neighbor_normalized_id": raw["qseqid"],
                "neighbor_strand": "+",
                "putative_gene_label": "aioB",
                "annotation": annotation,
                "qcov_pct": f"{100.0 * length / qlen:.2f}",
                "scov_pct": f"{100.0 * length / slen:.2f}",
                "ref_family": ref_family(raw["sseqid"]),
            }
            row["accepted_joint_reference_hit"] = "yes" if accepted_hit(row) else "no"
            if row["accepted_joint_reference_hit"] == "yes":
                output.append(row)
    return output


def normalize_hmm_scores() -> tuple[list[dict[str, str]], list[str]]:
    rows, fields = read_tsv(SOURCE_48 / "HMM_scores_49_source.tsv")
    for row in rows:
        protein_id = row["protein_id"]
        if "|" in protein_id and protein_id.split("|", 1)[0] in {"IS", "AS", "ES", "NS"}:
            row["protein_id"] = protein_id.split("|", 1)[1]
    rows.sort(key=lambda row: row["protein_id"])
    return rows, fields


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    gff = parse_gff(SOURCE_NEW / "R2111_N500_0-10_bin13.gff")
    kofam = parse_kofam(SOURCE_NEW / "R2111_N500_0-10_bin13.kofamscan.txt")
    if len(gff) != 4:
        raise SystemExit(f"Expected 4 genes on {CONTIG_ID}, found {len(gff)}")

    annotations = {gene_id: best_annotation(gene_id, kofam) for gene_id in gff}
    target_gene_index = 3
    new_long = []
    for gene_id, coord in sorted(gff.items(), key=lambda item: int(item[1]["gene_index"])):
        gene_index = int(coord["gene_index"])
        relative = gene_index - target_gene_index
        role = "target" if relative == 0 else ("upstream" if relative < 0 else "downstream")
        annotation, label = annotations[gene_id]
        new_long.append({
            "target_id": TARGET,
            "bin_id": BIN_ID,
            "contig_id": CONTIG_ID,
            "target_gene_index": str(target_gene_index),
            "neighbor_gene_index": str(gene_index),
            "relative_position": str(relative),
            "role": role,
            "neighbor_raw_id": gene_id,
            "neighbor_normalized_id": gene_id,
            "putative_gene_label": label,
            "annotation": annotation,
        })

    old_long, long_fields = read_tsv(SOURCE_48 / "gene_neighborhood_48_long.tsv")
    if TARGET in {row["target_id"] for row in old_long}:
        raise SystemExit("New target is already present in the 48-MAG long table")
    long_49 = old_long + new_long
    write_tsv(OUTPUT / "gene_neighborhood_49_long.tsv", long_49, long_fields)

    target_annotation = annotations[TARGET][0]
    neighborhood_notes = "; ".join(
        f"{row['relative_position']}:{row['annotation']}" for row in new_long
    )
    new_summary = {
        "target_id": TARGET,
        "bin_id": BIN_ID,
        "contig_id": CONTIG_ID,
        "target_gene_index": "3",
        "n_genes_in_window": "4",
        "target_annotation": target_annotation,
        "idrB_present": "No",
        "idrP1_present": "No",
        "idrP2_present": "No",
        "aioB_present": "Yes",
        "aio_operon_like": "Yes",
        "gene_neighborhood": "aioB/aio-like neighborhood",
        "neighborhood_notes": neighborhood_notes,
    }
    old_summary, summary_fields = read_tsv(SOURCE_48 / "gene_neighborhood_48_summary.tsv")
    summary_49 = old_summary + [new_summary]
    write_tsv(OUTPUT / "gene_neighborhood_49_summary.tsv", summary_49, summary_fields)

    hmm_49, hmm_fields = normalize_hmm_scores()
    write_tsv(OUTPUT / "HMM_scores_49.tsv", hmm_49, hmm_fields)
    hmm_by_id = {row["protein_id"]: row for row in hmm_49}
    hmm = hmm_by_id[TARGET]

    new_hits = parse_new_diamond_hits(annotations[f"{CONTIG_ID}_2"][0])
    old_hits, hit_fields = read_tsv(SOURCE_48 / "joint_reference_accepted_hits_long.tsv")
    hits_49 = old_hits + new_hits
    write_tsv(OUTPUT / "joint_reference_accepted_hits_long_49.tsv", hits_49, hit_fields)
    write_tsv(OUTPUT / "new_MAG_joint_reference_accepted_hits.tsv", new_hits, hit_fields)

    canonical_hits = [row for row in new_hits if row["ref_family"] == "canonical_AioB"]
    if not canonical_hits:
        raise SystemExit("Expected an accepted canonical AioB hit for the new MAG")
    best_b = max(canonical_hits, key=lambda row: float(row["bitscore"]))
    strict_fail = (
        "fewer_than_two_distinct_P_like_genes;overlapping_support_gene_coordinates;"
        "P_like_length_outside_150_700aa"
    )

    old_qc, qc_fields = read_tsv(SOURCE_48 / "strict_synteny_qc_48_reclassified.tsv")
    new_qc = {field: "" for field in qc_fields}
    new_qc.update({
        "candidate_id": TARGET,
        "contig_id": CONTIG_ID,
        "A_gene_id": TARGET,
        "A_start": gff[TARGET]["start"],
        "A_end": gff[TARGET]["end"],
        "A_strand": gff[TARGET]["strand"],
        "A_partial": gff[TARGET]["partial"],
        "B_hit_id": f"{CONTIG_ID}_2",
        "B_start": gff[f"{CONTIG_ID}_2"]["start"],
        "B_end": gff[f"{CONTIG_ID}_2"]["end"],
        "B_strand": gff[f"{CONTIG_ID}_2"]["strand"],
        "B_partial": gff[f"{CONTIG_ID}_2"]["partial"],
        "B_ref_family": "canonical_AioB",
        "B_reference_hit": best_b["sseqid"],
        "B_bitscore": best_b["bitscore"],
        "B_qcov_pct": best_b["qcov_pct"],
        "relative_gene_order": "B2(+)-A3(+)",
        "strand_pattern": "B:+;A:+",
        "distinct_P_like_count": "0",
        "same_contig": "yes",
        "maximum_gene_gap": "1",
        "arrow_plot_visible": "yes",
        "coordinates_non_overlapping": "no",
        "protein_lengths_reasonable": "no",
        "QC_pass": "no",
        "QC_fail_reasons": strict_fail,
        "strict_neighborhood_class": "partial DIRM-like neighborhood",
        "strict_neighborhood_class_revised": "partial DIRM-like neighborhood",
    })
    qc_49 = old_qc + [new_qc]
    write_tsv(OUTPUT / "strict_synteny_qc_49_reclassified.tsv", qc_49, qc_fields)

    old_final, final_fields = read_tsv(
        SOURCE_48 / "updated_final_evidence_table_joint_reference_strict_QC.tsv"
    )
    new_final = {field: "" for field in final_fields}
    new_final.update({
        "Sequence_ID": TARGET,
        "MAG_ID": BIN_ID,
        "Habitat": "NS",
        "Sample_group": "R2111_N500",
        "Combined_score": hmm["combined_bitscore"],
        "IdrA_score": hmm["idrA_bitscore"],
        "AioA_score": hmm["aioA_bitscore"],
        "IdrA_minus_AioA": hmm["idrA_minus_aioA"],
        "Combined_640_pass": "Yes" if float(hmm["combined_bitscore"]) >= 640 else "No",
        "HMM_direction": "IdrA-HMM higher",
        "Tree_clade": "Unknown-clade-associated",
        "Bootstrap_support": "100/100",
        "nearest_reference": NEAREST_REFERENCE,
        "nearest_reference_family": "Unknown_clade",
        "tree_distance": "1.138459",
        "gene_neighborhood": "aioB/aio-like neighborhood",
        "idrB_present": "No",
        "idrP1_present": "No",
        "idrP2_present": "No",
        "aioB_present": "Yes",
        "Final_class": "unknown/uncertain DMSOR",
        "Confidence": "Medium",
        "Evidence_summary": (
            f"IdrA-HMM higher; tree=Unknown-clade-associated; sister_MAG={SISTER}; "
            f"sister_node_support=100/100; nearest_ref={NEAREST_REFERENCE}; "
            "neighborhood=aioB/aio-like neighborhood"
        ),
        "Conflict_flags": (
            "IdrA HMM exceeds AioA HMM, but the canonical AioB neighborhood and unknown tree "
            "clade do not support definitive IdrA/DIRM assignment"
        ),
        "target_annotation": target_annotation,
        "neighborhood_notes": neighborhood_notes,
        "Previous_Final_class": "unknown/uncertain DMSOR",
        "Previous_Confidence": "Medium",
        "B_related_hit": f"{CONTIG_ID}_2",
        "B_related_similarity_or_tree_clade": "canonical_AioB",
        "P_like_count": "0",
        "P_like_hit_ids": "",
        "gene_order": "B2(+)-A3(+)",
        "strand_consistency": "B:+;A:+",
        "neighborhood_completeness": strict_fail,
        "neighborhood_class": "partial DIRM-like neighborhood",
        "synteny_supported_DIRM": "no",
        "strict_synteny_QC_pass": "no",
        "strict_synteny_QC_fail_reasons": strict_fail,
        "A_gene_coordinates": "1142-3616(+)",
        "B_gene_coordinates": "552-1142(+)",
        "P_like_1_coordinates": "-()",
        "P_like_2_coordinates": "-()",
        "classification_update": "retained_unknown_tree_boundary",
        "Joint_reference_evidence_summary": (
            "Strict synteny QC did not support a complete A-B-related-P-like-P-like "
            f"neighborhood (partial DIRM-like neighborhood); reasons: {strict_fail}."
        ),
    })
    final_49 = old_final + [new_final]
    write_tsv(OUTPUT / "updated_final_evidence_table_49_strict_QC.tsv", final_49, final_fields)

    checks = {
        "old_final_candidates": len(old_final),
        "final_49_candidates": len(final_49),
        "final_49_unique_ids": len({row["Sequence_ID"] for row in final_49}),
        "neighborhood_summary_candidates": len(summary_49),
        "neighborhood_summary_unique_ids": len({row["target_id"] for row in summary_49}),
        "neighborhood_long_rows": len(long_49),
        "neighborhood_long_unique_targets": len({row["target_id"] for row in long_49}),
        "strict_qc_candidates": len(qc_49),
        "hmm_score_candidates": len(hmm_49),
        "hmm_score_unique_ids": len({row["protein_id"] for row in hmm_49}),
        "new_target_accepted_joint_reference_hits": len(new_hits),
    }
    expected = {
        "old_final_candidates": 48,
        "final_49_candidates": 49,
        "final_49_unique_ids": 49,
        "neighborhood_summary_candidates": 49,
        "neighborhood_summary_unique_ids": 49,
        "neighborhood_long_unique_targets": 49,
        "strict_qc_candidates": 49,
        "hmm_score_candidates": 49,
        "hmm_score_unique_ids": 49,
    }
    for key, value in expected.items():
        if checks[key] != value:
            raise SystemExit(f"Validation failed: {key}={checks[key]}, expected {value}")
    with (OUTPUT / "validation_summary.txt").open("w") as handle:
        for key, value in checks.items():
            handle.write(f"{key}\t{value}\n")
        handle.write(f"new_target\t{TARGET}\n")
        handle.write(f"new_target_tree_label\t{TREE_LABEL}\n")
        handle.write(f"new_target_sister\t{SISTER}\n")
        handle.write("new_target_sister_support\t100/100\n")
        handle.write(f"new_target_nearest_reference\t{NEAREST_REFERENCE}\n")
        handle.write("new_target_nearest_reference_distance\t1.138459\n")


if __name__ == "__main__":
    main()
