#!/usr/bin/env python3
"""Integrate tree/synteny review and build auditable personalized references."""

import argparse
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path


CODONS = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    fields = fields or (list(rows[0]) if rows else [])
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def fasta(path):
    records, name, sequence = {}, None, []
    for line in Path(path).open(errors="replace"):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None: records[name] = "".join(sequence)
            name, sequence = line[1:].split()[0], []
        elif name is not None: sequence.append(line)
    if name is not None: records[name] = "".join(sequence)
    return records


def write_fasta(path, records):
    with Path(path).open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def translate(sequence):
    sequence = sequence.upper().replace("U", "T")
    return "".join(CODONS.get(sequence[i:i + 3], "X") for i in range(0, len(sequence) - 2, 3))


def translation_matches(cds_sequence, protein_sequence):
    observed = translate(cds_sequence).rstrip("*")
    expected = protein_sequence.rstrip("*")
    if expected.startswith("M") and observed:
        observed = "M" + observed[1:]
    return observed == expected and "*" not in observed


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""): h.update(chunk)
    return h.hexdigest()


def bin_index(bin_dir):
    assignments, bin_paths, bin_records = defaultdict(list), {}, {}
    for path in sorted(Path(bin_dir).glob("*")):
        if not path.is_file() or path.suffix.lower() not in {".fa", ".fna", ".fasta"}: continue
        mag = path.stem; records = fasta(path)
        bin_paths[mag], bin_records[mag] = str(path.resolve()), records
        for contig in records: assignments[contig].append(mag)
    return assignments, bin_paths, bin_records


def classify(phylogeny, complete_synteny):
    if phylogeny == "phylogenetically strict-core IdrA-associated":
        return "strict synteny-supported DIRM-like IdrA-associated" if complete_synteny else "partial IdrA-associated"
    if phylogeny == "phylogenetically partial IdrA-associated": return "partial IdrA-associated"
    if phylogeny == "canonical AioA": return "canonical AioA"
    if phylogeny == "other DMSOR-family": return "other DMSOR-family"
    return "unresolved/truncated"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--review", required=True)
    parser.add_argument("--neighborhood", required=True)
    parser.add_argument("--faa", required=True)
    parser.add_argument("--fna", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--bin-dir", required=True)
    parser.add_argument("--mag-metadata", help="Optional TSV with an explicit MAG_id column")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args(); out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    review, neighborhood = read_tsv(args.review), read_tsv(args.neighborhood)
    proteins, cds, contigs = fasta(args.faa), fasta(args.fna), fasta(args.contigs)
    assignments, bin_paths, bin_records = bin_index(args.bin_dir)
    allowed_phylogeny = {
        "phylogenetically strict-core IdrA-associated", "phylogenetically partial IdrA-associated",
        "canonical AioA", "other DMSOR-family", "unresolved/truncated",
    }
    if any(row["final_phylogenetic_class"] not in allowed_phylogeny for row in review):
        raise SystemExit("Tree review is incomplete or contains a non-permitted phylogenetic class")
    neighbors_by_candidate = defaultdict(list)
    for row in neighborhood: neighbors_by_candidate[row["candidate_id"]].append(row)

    final_rows, gene_rows, cluster_rows, candidate_mag_rows = [], [], [], []
    individual_cds, individual_faa, cluster_fasta = [], [], []
    candidate_contig_ids = defaultdict(list)
    accepted_strict = []
    for candidate in review:
        candidate_id, contig_id = candidate["protein_id"], candidate["contig_id"]
        genes = neighbors_by_candidate[candidate_id]
        focal = next((row for row in genes if row["neighbor_id"] == candidate_id), None)
        if not focal: raise SystemExit(f"Missing focal neighborhood row: {candidate_id}")
        strand = focal["strand"]; direction = 1 if strand == "+" else -1
        accepted = [row for row in genes if row["neighbor_annotation_accepted"].lower() == "yes"]
        b_hits = [row for row in accepted if row["neighbor_family"] == "IdrB_related"]
        p_hits = [row for row in accepted if row["neighbor_family"] == "P_like"]
        b_hits.sort(key=lambda row: abs(int(row["relative_ORF"])))
        p_hits.sort(key=lambda row: direction * int(row["relative_ORF"]))
        roles = [("IdrA", focal)]
        if b_hits: roles.append(("IdrB", b_hits[0]))
        downstream_p = [row for row in p_hits if direction * int(row["relative_ORF"]) > 0]
        if downstream_p: roles.append(("P-like-1", downstream_p[0]))
        if len(downstream_p) > 1: roles.append(("P-like-2", downstream_p[1]))
        role_positions = {role: direction * int(row["relative_ORF"]) for role, row in roles}
        compact_order = all(role_positions.get(role) == position for position, role in enumerate(("IdrA","IdrB","P-like-1","P-like-2")))
        same_strand = len({row["strand"] for _, row in roles}) == 1
        nonoverlap = all(int(left[1]["end"]) < int(right[1]["start"]) for left, right in zip(sorted(roles, key=lambda x:int(x[1]["start"])), sorted(roles, key=lambda x:int(x[1]["start"]))[1:])) if len(roles) == 4 else False
        complete_synteny = len(roles) == 4 and compact_order and same_strand and nonoverlap and focal["edge_censored"] == "no"
        final_class = classify(candidate["final_phylogenetic_class"], complete_synteny)
        mags = assignments.get(contig_id, [])
        mag_status = "unbinned" if not mags else "unique" if len(mags) == 1 else "ambiguous_multiple_bins"
        mag = mags[0] if len(mags) == 1 else ""
        final = dict(candidate)
        final.update({"IdrB_present": "yes" if b_hits else "no", "P_like_count": len(p_hits),
                      "compact_A_B_P-like-1_P-like-2_order": "yes" if compact_order else "no",
                      "same_strand": "yes" if same_strand else "no", "nonoverlapping": "yes" if nonoverlap else "no",
                      "edge_censored": focal["edge_censored"], "final_integrated_class": final_class,
                      "MAG_assignment_status": mag_status, "MAG_id": mag})
        final_rows.append(final)
        candidate_mag_rows.append({"candidate_id": candidate_id, "contig_id": contig_id,
                                   "MAG_assignment_status": mag_status, "MAG_id": mag,
                                   "all_matching_MAGs": ";".join(mags),
                                   "final_integrated_class": final_class})
        if final_class != "strict synteny-supported DIRM-like IdrA-associated": continue
        accepted_strict.append(candidate_id)
        role_rows = []
        for role, row in roles:
            gene_id = row["neighbor_id"]
            if gene_id not in cds or gene_id not in proteins:
                raise SystemExit(f"Missing CDS/protein for strict-cluster gene {gene_id}")
            translation_status = "PASS_exact" if translation_matches(cds[gene_id], proteins[gene_id]) else "FAIL"
            if translation_status != "PASS_exact":
                raise SystemExit(f"Translation validation failed for {gene_id}")
            header = f"22_N30_16|{candidate_id}|{role}|contig={contig_id}|protein={gene_id}|MAG={mag or 'unbinned'}"
            individual_cds.append((header, cds[gene_id])); individual_faa.append((header, proteins[gene_id]))
            gene = {"candidate_id": candidate_id, "role": role, "gene_id": gene_id,
                    "contig_id": contig_id, "start": row["start"], "end": row["end"],
                    "strand": row["strand"], "MAG_id": mag, "translation_status": translation_status,
                    "reference_id": header}
            gene_rows.append(gene); role_rows.append(gene)
        genomic = sorted(role_rows, key=lambda row: int(row["start"]))
        region_start, region_end = int(genomic[0]["start"]), int(genomic[-1]["end"])
        region_id = f"22_N30_16|{candidate_id}|A_B_P-like-1_P-like-2_native_cluster|contig={contig_id}|MAG={mag or 'unbinned'}"
        cluster_fasta.append((region_id, contigs[contig_id][region_start-1:region_end]))
        candidate_contig_ids[(contig_id, mag or "unbinned")].append(candidate_id)
        role_lookup = {row["role"]: row for row in role_rows}
        for first_role, second_role in (("IdrA","IdrB"),("IdrB","P-like-1"),("P-like-1","P-like-2")):
            first, second = role_lookup[first_role], role_lookup[second_role]
            genomic_pair = sorted((first, second), key=lambda row: int(row["start"]))
            left_gene, right_gene = genomic_pair
            cluster_rows.append({"candidate_id": candidate_id, "contig_id": contig_id,
                "first_role_transcriptional": first_role, "second_role_transcriptional": second_role,
                "first_gene_start": first["start"], "first_gene_end": first["end"],
                "second_gene_start": second["start"], "second_gene_end": second["end"],
                "genomic_left_gene_end": left_gene["end"], "genomic_right_gene_start": right_gene["start"],
                "intergenic_bp": int(right_gene["start"])-int(left_gene["end"])-1,
                "strand": strand, "MAG_id": mag})

    write_tsv(out/"22_N30_16_final_candidate_classification.tsv", final_rows)
    write_tsv(out/"22_N30_16_candidate_to_MAG.tsv", candidate_mag_rows)
    metadata = {}
    metadata_fields = []
    if args.mag_metadata:
        metadata_rows = read_tsv(args.mag_metadata)
        if metadata_rows and "MAG_id" not in metadata_rows[0]:
            raise SystemExit("--mag-metadata must contain an exact MAG_id column")
        metadata = {row["MAG_id"]: row for row in metadata_rows}
        metadata_fields = [field for field in (list(metadata_rows[0]) if metadata_rows else []) if field != "MAG_id"]
    by_mag = defaultdict(list)
    for row in candidate_mag_rows:
        if row["MAG_id"]: by_mag[row["MAG_id"]].append(row)
    host_rows = []
    for mag, candidates in sorted(by_mag.items()):
        host_rows.append({"MAG_id": mag, "MAG_FASTA": bin_paths.get(mag, ""),
                          "candidate_count": len(candidates),
                          "candidate_ids": ";".join(row["candidate_id"] for row in candidates),
                          "candidate_classes": ";".join(row["final_integrated_class"] for row in candidates),
                          **metadata.get(mag, {})})
    write_tsv(out/"22_N30_16_candidate_host_MAG_summary.tsv", host_rows,
              ["MAG_id","MAG_FASTA","candidate_count","candidate_ids","candidate_classes"] + metadata_fields)
    write_tsv(out/"22_N30_16_strict_cluster_genes.tsv", gene_rows,
              ["candidate_id","role","gene_id","contig_id","start","end","strand","MAG_id","translation_status","reference_id"])
    write_tsv(out/"22_N30_16_strict_cluster_boundaries.tsv", cluster_rows,
              ["candidate_id","contig_id","first_role_transcriptional","second_role_transcriptional",
               "first_gene_start","first_gene_end","second_gene_start","second_gene_end",
               "genomic_left_gene_end","genomic_right_gene_start","intergenic_bp","strand","MAG_id"])
    write_fasta(out/"22_N30_16_strict_cluster_individual_CDS.fna", individual_cds)
    write_fasta(out/"22_N30_16_strict_cluster_individual_proteins.faa", individual_faa)
    write_fasta(out/"22_N30_16_strict_A_B_P-like-1_P-like-2_native_regions.fna", cluster_fasta)
    write_fasta(out/"22_N30_16_strict_candidate_contigs.fna", [
        (f"22_N30_16|candidates={','.join(ids)}|full_contig|contig={contig_id}|MAG={mag}", contigs[contig_id])
        for (contig_id, mag), ids in sorted(candidate_contig_ids.items())
    ])
    write_fasta(out/"22_N30_16_all_K08356_candidates_CDS.fna", [
        (f"22_N30_16|{row['protein_id']}|{row['final_integrated_class'].replace(' ', '_')}", cds[row["protein_id"]])
        for row in final_rows
    ])

    all_bins, bin_manifest = [], []
    host_mags = {row["MAG_id"] for row in candidate_mag_rows if row["MAG_id"]}
    host_bins = []
    for mag, records in sorted(bin_records.items()):
        for contig, sequence in records.items():
            ref_id = f"MAG={mag}|contig={contig}"
            all_bins.append((ref_id, sequence)); bin_manifest.append({"reference_id": ref_id, "MAG_id": mag, "original_contig_id": contig, "length": len(sequence)})
            if mag in host_mags: host_bins.append((ref_id, sequence))
    write_fasta(out/"22_N30_16_all_bins_competitive.fna", all_bins)
    write_fasta(out/"22_N30_16_candidate_host_MAGs.fna", host_bins)
    write_tsv(out/"22_N30_16_all_bins_reference_manifest.tsv", bin_manifest,
              ["reference_id","MAG_id","original_contig_id","length"])

    generated = sorted(path for path in out.iterdir() if path.is_file())
    with (out/"checksums.sha256").open("w") as handle:
        for path in generated: handle.write(f"{sha256(path)}  {path.name}\n")
    (out/"FINALIZATION_SUMMARY.txt").write_text(
        f"T640 candidates\t{len(final_rows)}\nstrict synteny-supported clusters\t{len(accepted_strict)}\n"
        "Four CDS are retained as separate references; they are not concatenated for gene-level quantification.\n"
    )


if __name__ == "__main__": main()
