#!/usr/bin/env python3
"""Integrate tree, neighborhood, bin, and MAG evidence for 82 public candidates."""

import argparse
import csv
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo


SAMPLES = ("22_N10_3", "22_N30_16", "23_N10_9")
EXPECTED_COUNTS = {"22_N10_3": 22, "22_N30_16": 10, "23_N10_9": 50}
ANCHOR_CLASS = {
    "IdrA clade": "strict_core_IdrA_associated",
    "IdrA-associated clade": "partial_IdrA_associated",
    "Canonical AioA clade": "canonical_AioA",
    "Unknown DMSOR clade": "other_DMSOR",
}


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    fields = fields or (list(rows[0]) if rows else [])
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def iter_fasta(path):
    name = description = None
    sequence = []
    with Path(path).open(errors="replace") as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, description, "".join(sequence)
                description = line[1:]
                name = description.split()[0]
                sequence = []
            elif name is not None:
                sequence.append(line)
    if name is not None:
        yield name, description, "".join(sequence)


def write_fasta(path, records):
    with Path(path).open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def parse_tbl(path):
    hits = {}
    with Path(path).open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            hits[fields[0]] = {"evalue": fields[4], "bitscore": float(fields[5])}
    return hits


def parse_domtbl(path):
    hits = defaultdict(list)
    with Path(path).open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            target_length = int(fields[2])
            start, end = int(fields[17]), int(fields[18])
            coverage = (end - start + 1) / target_length if target_length else 0
            hits[fields[0]].append((coverage, start, end))
    return hits


def parse_prodigal_header(sequence_id, description):
    pieces = [piece.strip() for piece in description.split("#")]
    if len(pieces) < 4:
        return None
    try:
        start, end, strand_code = int(pieces[1]), int(pieces[2]), int(pieces[3])
    except ValueError:
        return None
    return sequence_id.rsplit("_", 1)[0], start, end, "+" if strand_code == 1 else "-"


def parse_gff_targets(path, target_contigs):
    genes = defaultdict(list)
    with Path(path).open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS" or fields[0] not in target_contigs:
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            genes[fields[0]].append({
                "gff_id": attrs.get("ID", ""), "start": int(fields[3]),
                "end": int(fields[4]), "strand": fields[6],
            })
    for rows in genes.values():
        rows.sort(key=lambda row: (row["start"], row["end"]))
    return genes


def load_target_gene_sequences(path, target_contigs):
    records, descriptions, coordinates = {}, {}, {}
    for sequence_id, description, sequence in iter_fasta(path):
        if sequence_id.rsplit("_", 1)[0] not in target_contigs:
            continue
        records[sequence_id] = sequence
        descriptions[sequence_id] = description
        coordinate = parse_prodigal_header(sequence_id, description)
        if coordinate:
            coordinates[coordinate] = sequence_id
    return records, descriptions, coordinates


def load_target_contigs(path, target_contigs):
    return {name: sequence for name, _, sequence in iter_fasta(path) if name in target_contigs}


def sample_inputs(pilot_root, sample):
    sample_root = pilot_root / sample
    hmm = sample_root / "04_hmm"
    prodigal = sample_root / "03_prodigal"
    return {
        "combined": hmm / "combined_iriA_aioA.raw.tblout",
        "strict": hmm / "combined_iriA_aioA.strict_T640.tblout",
        "domtbl": hmm / "combined_iriA_aioA.raw.domtblout",
        "iria": hmm / "iriA_new.raw.tblout",
        "aioa": hmm / "aioA.raw.tblout",
        "faa": prodigal / f"{sample}.proteins.faa",
        "fna": prodigal / f"{sample}.genes.fna",
        "gff": prodigal / f"{sample}.genes.gff",
        "contigs": sample_root / "02_megahit" / "final.contigs.fa",
    }


def build_candidate_inputs(pilot_root, outdir):
    candidate_rows, neighborhoods, neighborhood_proteins = [], [], []
    for sample in SAMPLES:
        paths = sample_inputs(pilot_root, sample)
        for label, path in paths.items():
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f"Missing {sample} {label}: {path}")
        combined = parse_tbl(paths["combined"])
        strict = parse_tbl(paths["strict"])
        iria = parse_tbl(paths["iria"])
        aioa = parse_tbl(paths["aioa"])
        domains = parse_domtbl(paths["domtbl"])
        if len(strict) != EXPECTED_COUNTS[sample]:
            raise SystemExit(f"{sample}: expected {EXPECTED_COUNTS[sample]} strict candidates, found {len(strict)}")
        target_contigs = {candidate.rsplit("_", 1)[0] for candidate in strict}
        genes = parse_gff_targets(paths["gff"], target_contigs)
        proteins, descriptions, coordinate_ids = load_target_gene_sequences(paths["faa"], target_contigs)
        cds, _, _ = load_target_gene_sequences(paths["fna"], target_contigs)
        contigs = load_target_contigs(paths["contigs"], target_contigs)
        for contig_id, rows in genes.items():
            for row in rows:
                row["protein_id"] = coordinate_ids.get((contig_id, row["start"], row["end"], row["strand"]), "")
        for candidate in sorted(strict, key=lambda item: strict[item]["bitscore"], reverse=True):
            contig_id = candidate.rsplit("_", 1)[0]
            if candidate not in proteins or candidate not in cds or contig_id not in contigs:
                raise SystemExit(f"Failed exact sequence extraction for {sample}|{candidate}")
            gene_list = genes[contig_id]
            index = next((i for i, row in enumerate(gene_list) if row["protein_id"] == candidate), None)
            if index is None:
                raise SystemExit(f"Failed exact GFF reconciliation for {sample}|{candidate}")
            gene = gene_list[index]
            contig_length = len(contigs[contig_id])
            edge_distance = min(gene["start"] - 1, contig_length - gene["end"])
            best_domain = max(domains.get(candidate, [(0, "", "")]), key=lambda item: item[0])
            description = descriptions[candidate]
            partial_match = re.search(r"partial=([01]{2})", description)
            candidate_rows.append({
                "sample": sample, "candidate_id": f"{sample}|{candidate}", "protein_id": candidate,
                "contig_id": contig_id, "gene_start": gene["start"], "gene_end": gene["end"],
                "strand": gene["strand"], "protein_length_aa": len(proteins[candidate].rstrip("*")),
                "CDS_length_nt": len(cds[candidate]), "contig_length_bp": contig_length,
                "distance_to_contig_end_bp": edge_distance,
                "prodigal_partial_code": partial_match.group(1) if partial_match else "",
                "combined_bitscore": combined[candidate]["bitscore"],
                "combined_evalue": combined[candidate]["evalue"],
                "iriA_bitscore": iria.get(candidate, {}).get("bitscore", ""),
                "iriA_evalue": iria.get(candidate, {}).get("evalue", ""),
                "aioA_bitscore": aioa.get(candidate, {}).get("bitscore", ""),
                "aioA_evalue": aioa.get(candidate, {}).get("evalue", ""),
                "combined_HMM_coverage": round(best_domain[0], 6),
                "combined_HMM_alignment_start_aa": best_domain[1],
                "combined_HMM_alignment_end_aa": best_domain[2],
            })
            edge_censored = index < 10 or len(gene_list) - index - 1 < 10
            for neighbor_index in range(max(0, index - 10), min(len(gene_list), index + 11)):
                neighbor = gene_list[neighbor_index]
                neighbor_id = neighbor["protein_id"]
                query_id = f"{sample}|{neighbor_id}" if neighbor_id else ""
                neighborhoods.append({
                    "sample": sample, "candidate_id": f"{sample}|{candidate}",
                    "protein_id": candidate, "relative_ORF": neighbor_index - index,
                    "contig_id": contig_id, "neighbor_id": neighbor_id,
                    "query_id": query_id, "start": neighbor["start"], "end": neighbor["end"],
                    "strand": neighbor["strand"], "contig_length_bp": contig_length,
                    "edge_censored": "yes" if edge_censored else "no",
                })
                if neighbor_id and neighbor_id in proteins:
                    neighborhood_proteins.append((query_id, proteins[neighbor_id]))
    if len(candidate_rows) != 82:
        raise SystemExit(f"Expected 82 candidates, found {len(candidate_rows)}")
    write_tsv(outdir / "external82_candidate_base_evidence.tsv", candidate_rows)
    write_tsv(outdir / "external82_neighborhood_raw.tsv", neighborhoods)
    deduplicated = dict(neighborhood_proteins)
    write_fasta(outdir / "external82_neighborhood_proteins.faa", sorted(deduplicated.items()))
    return candidate_rows, neighborhoods


def concatenate_refs(idrb_ref, plike_ref, output):
    records = []
    for path in (idrb_ref, plike_ref):
        records.extend((name, sequence) for name, _, sequence in iter_fasta(path))
    write_fasta(output, records)


def run_diamond(diamond, query, refs, outdir, threads):
    combined_refs = outdir / "validated_IdrB_P_like_refs.faa"
    concatenate_refs(refs[0], refs[1], combined_refs)
    database = outdir / "validated_IdrB_P_like_refs"
    hits = outdir / "external82_validated_neighbor_hits.raw.tsv"
    subprocess.run([diamond, "makedb", "--in", str(combined_refs), "-d", str(database)], check=True,
                   stdout=subprocess.DEVNULL)
    subprocess.run([
        diamond, "blastp", "--query", str(query), "--db", str(database) + ".dmnd",
        "--out", str(hits), "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qlen", "slen", "--more-sensitive", "--evalue", "1e-5", "--max-target-seqs", "25",
        "--threads", str(threads),
    ], check=True)
    return hits


def reference_family(reference):
    if reference.startswith(("IdrP_like", "IdrP1", "IdrP2")):
        return "P_like"
    if reference.startswith(("IdrB", "AioB")):
        return "IdrB_related"
    return "other"


def accepted_neighbor(group, evalue, bitscore, qcov, scov):
    if group == "IdrB_related":
        return evalue <= 1e-5 and bitscore >= 50 and qcov >= 50 and scov >= 35
    if group == "P_like":
        return evalue <= 1e-10 and bitscore >= 80 and qcov >= 45 and scov >= 35
    return False


def annotate_neighbors(neighborhoods, hits_path, outdir):
    best, hit_rows = {}, []
    with Path(hits_path).open(newline="") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if len(fields) < 14:
                continue
            query, reference = fields[0], fields[1]
            alignment_length = float(fields[3])
            evalue, bitscore = float(fields[10]), float(fields[11])
            qlen, slen = float(fields[12]), float(fields[13])
            group = reference_family(reference)
            qcov = 100 * alignment_length / qlen if qlen else 0
            scov = 100 * alignment_length / slen if slen else 0
            keep = accepted_neighbor(group, evalue, bitscore, qcov, scov)
            row = {
                "query_id": query, "reference_id": reference, "family": group,
                "evalue": evalue, "bitscore": bitscore,
                "query_coverage_pct": round(qcov, 3), "subject_coverage_pct": round(scov, 3),
                "accepted": "yes" if keep else "no",
            }
            hit_rows.append(row)
            if keep and (query not in best or bitscore > best[query]["bitscore"]):
                best[query] = row
    reviewed = []
    for row in neighborhoods:
        hit = best.get(row["query_id"])
        reviewed.append({
            **row,
            "neighbor_family": hit["family"] if hit else "",
            "neighbor_annotation_accepted": "yes" if hit else "no",
            "neighbor_bitscore": hit["bitscore"] if hit else "",
            "neighbor_reference": hit["reference_id"] if hit else "",
        })
    write_tsv(outdir / "external82_neighborhood_reviewed.tsv", reviewed)
    write_tsv(outdir / "external82_validated_neighbor_hits.tsv", hit_rows)
    return reviewed


def normalize_anchor_protein(protein_id):
    return re.sub(r"_bin[0-9]+(?=-k141_)", "", protein_id)


def anchor_leaves(tree_fasta, evidence):
    q_headers = [name for name, _, _ in iter_fasta(tree_fasta) if name.startswith("Q")]
    by_suffix = {name.split("__", 1)[1]: name for name in q_headers if "__" in name}
    anchors, missing = {}, []
    for row in evidence:
        if row["source_class"] != "core49_T640_intersection":
            continue
        normalized = normalize_anchor_protein(row["protein_ID"])
        leaf = by_suffix.get(normalized)
        if not leaf:
            missing.append(row["protein_ID"])
            continue
        anchors[leaf] = {
            "protein_id": row["protein_ID"], "final_class": row["Final_class"],
            "anchor_class": ANCHOR_CLASS[row["Final_class"]],
            "strict_synteny": row["strict_synteny_QC_pass"],
        }
    if len(anchors) != 48 or missing:
        raise SystemExit(f"Expected 48 mapped reviewed core anchors; found {len(anchors)}; missing={missing}")
    return anchors


def tree_evidence(tree_path, tree_fasta, anchor_evidence, candidate_rows):
    tree = Phylo.read(tree_path, "newick")
    terminals = {terminal.name: terminal for terminal in tree.get_terminals()}
    anchors = anchor_leaves(tree_fasta, read_tsv(anchor_evidence))
    anchor_nodes = [(name, terminals[name], info) for name, info in anchors.items()]
    result = {}
    for row in candidate_rows:
        candidate_id = row["candidate_id"]
        if candidate_id not in terminals:
            raise SystemExit(f"Candidate absent from tree: {candidate_id}")
        node = terminals[candidate_id]
        ranked = sorted((tree.distance(node, anchor_node), name, anchor_node, info)
                        for name, anchor_node, info in anchor_nodes)
        nearest_distance, nearest_name, nearest_node, nearest_info = ranked[0]
        class_minimum = {}
        for distance, _, _, info in ranked:
            class_minimum.setdefault(info["anchor_class"], distance)
        class_rank = sorted(class_minimum.items(), key=lambda item: item[1])
        second_class = class_rank[1] if len(class_rank) > 1 else ("", float("nan"))
        ancestor = tree.common_ancestor(node, nearest_node)
        descendant_names = {terminal.name for terminal in ancestor.get_terminals()}
        descendant_classes = sorted({anchors[name]["anchor_class"] for name in descendant_names if name in anchors})
        support = ancestor.confidence if ancestor.confidence is not None else ""
        monophyletic = len(descendant_classes) == 1
        provisional = nearest_info["anchor_class"] if monophyletic else "unresolved_tree_placement"
        margin = second_class[1] - nearest_distance
        review_reasons = []
        if not monophyletic:
            review_reasons.append("mixed_anchor_classes_below_MRCA")
        if support != "" and float(support) < 70:
            review_reasons.append("MRCA_support_below_70")
        if margin < 0.05:
            review_reasons.append("nearest_class_margin_below_0.05")
        manual = "yes" if review_reasons else "no"
        result[candidate_id] = {
            "nearest_anchor_leaf": nearest_name,
            "nearest_anchor_protein": nearest_info["protein_id"],
            "nearest_anchor_reviewed_class": nearest_info["final_class"],
            "nearest_anchor_distance": round(nearest_distance, 8),
            "second_nearest_anchor_class": second_class[0],
            "second_class_distance_margin": round(margin, 8),
            "candidate_anchor_MRCA_support": support,
            "anchor_classes_below_MRCA": ";".join(descendant_classes),
            "tree_monophyletic_anchor_class": "yes" if monophyletic else "no",
            "provisional_tree_class": provisional,
            "tree_manual_review_required": manual,
            "tree_review_reason": ";".join(review_reasons),
        }
    return result


def taxonomy_parts(classification):
    parts = {item.split("__", 1)[0]: item for item in classification.split(";") if "__" in item}
    return {
        "GTDB_phylum": parts.get("p", ""), "GTDB_class": parts.get("c", ""),
        "GTDB_order": parts.get("o", ""), "GTDB_family": parts.get("f", ""),
        "GTDB_genus": parts.get("g", ""),
    }


def integrate(candidate_rows, reviewed_neighbors, tree_rows, contig_bin_path, mag_summary_path, taxonomy_path):
    neighbors = defaultdict(list)
    for row in reviewed_neighbors:
        neighbors[row["candidate_id"]].append(row)
    contig_bins = {(row["sample"], row["protein_id"]): row for row in read_tsv(contig_bin_path)}
    mag_summary = {(row["sample"], row["bin_id"]): row for row in read_tsv(mag_summary_path)}
    taxonomy = {(row["sample"], row["bin_id"]): row for row in read_tsv(taxonomy_path)}
    final_rows = []
    for candidate in candidate_rows:
        candidate_id = candidate["candidate_id"]
        focal = next(row for row in neighbors[candidate_id] if int(row["relative_ORF"]) == 0)
        accepted = [row for row in neighbors[candidate_id] if row["neighbor_annotation_accepted"] == "yes"]
        positions = defaultdict(list)
        for row in accepted:
            positions[row["neighbor_family"]].append(int(row["relative_ORF"]))
        role_sets = []
        for side in (1, -1):
            b_exact = [row for row in accepted if row["neighbor_family"] == "IdrB_related" and int(row["relative_ORF"]) == side]
            p1_exact = [row for row in accepted if row["neighbor_family"] == "P_like" and int(row["relative_ORF"]) == 2 * side]
            p2_exact = [row for row in accepted if row["neighbor_family"] == "P_like" and int(row["relative_ORF"]) == 3 * side]
            role_sets.append((side, [focal] + b_exact[:1] + p1_exact[:1] + p2_exact[:1]))
        cluster_side, roles = next(((side, role_set) for side, role_set in role_sets if len(role_set) == 4), (0, [focal]))
        same_strand = len(roles) == 4 and len({row["strand"] for row in roles}) == 1
        genomic = sorted(roles, key=lambda row: int(row["start"]))
        nonoverlap = len(genomic) == 4 and all(int(left["end"]) < int(right["start"]) for left, right in zip(genomic, genomic[1:]))
        complete = len(roles) == 4 and same_strand and nonoverlap
        tree = tree_rows[candidate_id]
        phylogeny = tree["provisional_tree_class"]
        if tree["tree_manual_review_required"] == "yes":
            final_class = "unresolved_or_truncated"
        elif phylogeny == "strict_core_IdrA_associated" and complete:
            final_class = "strict_synteny-supported_DIRM-like_IdrA-associated"
        elif phylogeny in {"strict_core_IdrA_associated", "partial_IdrA_associated"}:
            final_class = "partial_IdrA-associated"
        elif phylogeny == "canonical_AioA":
            final_class = "canonical_AioA"
        elif phylogeny == "other_DMSOR":
            final_class = "other_DMSOR"
        else:
            final_class = "unresolved_or_truncated"
        bin_row = contig_bins.get((candidate["sample"], candidate["protein_id"]), {})
        bin_id = bin_row.get("bin_id", "") if bin_row.get("bin_status") == "binned" else ""
        mag = mag_summary.get((candidate["sample"], bin_id), {})
        tax = taxonomy.get((candidate["sample"], bin_id), {})
        classification = tax.get("gtdb_classification", "")
        final_rows.append({
            **candidate, **tree,
            "IdrB_related_relative_ORFs": ";".join(map(str, sorted(positions["IdrB_related"]))),
            "P_like_relative_ORFs": ";".join(map(str, sorted(positions["P_like"]))),
            "IdrB_related_present": "yes" if positions["IdrB_related"] else "no",
            "P_like_count": len(positions["P_like"]),
            "compact_A_B_P-like-1_P-like-2": "yes" if complete else "no",
            "cluster_relative_side": cluster_side if complete else "",
            "same_strand_complete_cluster": "yes" if same_strand else "no",
            "nonoverlapping_complete_cluster": "yes" if nonoverlap else "no",
            "edge_censored": focal["edge_censored"],
            "final_integrated_class": final_class,
            "final_manual_review_required": tree["tree_manual_review_required"],
            "bin_status": bin_row.get("bin_status", "unbinned"), "bin_id": bin_id,
            "MAG_completeness": mag.get("Completeness", ""),
            "MAG_contamination": mag.get("Contamination", ""),
            "MAG_quality_class": mag.get("quality_class", ""),
            "MAG_CoverM_mean": mag.get("CoverM_Mean", ""),
            "MAG_relative_abundance_pct": mag.get("CoverM_Relative_Abundance_pct", ""),
            "GTDB_classification": classification, **taxonomy_parts(classification),
            "GTDB_warnings": tax.get("gtdb_warnings", ""),
        })
    return final_rows


def write_summaries(rows, outdir):
    class_counts = Counter((row["sample"], row["final_integrated_class"]) for row in rows)
    classes = [
        "strict_synteny-supported_DIRM-like_IdrA-associated", "partial_IdrA-associated",
        "canonical_AioA", "other_DMSOR", "unresolved_or_truncated",
    ]
    summary = []
    for sample in SAMPLES:
        for final_class in classes:
            summary.append({"sample": sample, "final_integrated_class": final_class,
                            "candidate_count": class_counts[(sample, final_class)]})
    for final_class in classes:
        summary.append({"sample": "ALL", "final_integrated_class": final_class,
                        "candidate_count": sum(class_counts[(sample, final_class)] for sample in SAMPLES)})
    write_tsv(outdir / "external82_class_summary.tsv", summary)

    mag_groups = defaultdict(list)
    for row in rows:
        if row["bin_id"]:
            mag_groups[(row["sample"], row["bin_id"])].append(row)
    mag_rows = []
    for (sample, bin_id), candidates in sorted(mag_groups.items()):
        first = candidates[0]
        mag_rows.append({
            "sample": sample, "bin_id": bin_id, "candidate_count": len(candidates),
            "candidate_ids": ";".join(row["protein_id"] for row in candidates),
            "candidate_classes": ";".join(sorted({row["final_integrated_class"] for row in candidates})),
            "Completeness": first["MAG_completeness"], "Contamination": first["MAG_contamination"],
            "quality_class": first["MAG_quality_class"], "CoverM_Mean": first["MAG_CoverM_mean"],
            "relative_abundance_pct": first["MAG_relative_abundance_pct"],
            "GTDB_classification": first["GTDB_classification"], "GTDB_family": first["GTDB_family"],
            "GTDB_warnings": first["GTDB_warnings"],
        })
    write_tsv(outdir / "external82_candidate_host_MAG_summary.tsv", mag_rows)

    total = len(rows)
    lines = [
        "# External 82-candidate terminal classification", "",
        f"- Input candidates: {total} combined-HMM T640 K08356/DMSOR candidates.",
        "- These candidates were not pre-labelled as IdrA.",
        f"- Binned candidates: {sum(row['bin_status'] == 'binned' for row in rows)}.",
        f"- Candidates requiring tree review: {sum(row['final_manual_review_required'] == 'yes' for row in rows)}.",
        "", "## Integrated classes", "", "| Class | Count |", "|---|---:|",
    ]
    totals = Counter(row["final_integrated_class"] for row in rows)
    lines.extend(f"| {final_class} | {totals[final_class]} |" for final_class in classes)
    lines.extend(["", "Strict IdrA-associated status requires both strict-core tree placement and the complete validated DIRM-like neighborhood.", ""])
    (outdir / "EXTERNAL82_TERMINAL_CLASSIFICATION_REPORT.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-root", required=True, type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--tree-fasta", required=True, type=Path)
    parser.add_argument("--anchor-evidence", required=True, type=Path)
    parser.add_argument("--idrb-ref", required=True, type=Path)
    parser.add_argument("--plike-ref", required=True, type=Path)
    parser.add_argument("--contig-bin", required=True, type=Path)
    parser.add_argument("--mag-summary", required=True, type=Path)
    parser.add_argument("--mag-taxonomy", required=True, type=Path)
    parser.add_argument("--diamond", default="diamond")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    candidate_rows, neighborhoods = build_candidate_inputs(args.pilot_root, args.outdir)
    hits = run_diamond(args.diamond, args.outdir / "external82_neighborhood_proteins.faa",
                       (args.idrb_ref, args.plike_ref), args.outdir, args.threads)
    reviewed = annotate_neighbors(neighborhoods, hits, args.outdir)
    tree_rows = tree_evidence(args.tree, args.tree_fasta, args.anchor_evidence, candidate_rows)
    final = integrate(candidate_rows, reviewed, tree_rows, args.contig_bin, args.mag_summary, args.mag_taxonomy)
    write_tsv(args.outdir / "external82_integrated_candidate_classification.tsv", final)
    write_summaries(final, args.outdir)


if __name__ == "__main__":
    main()
