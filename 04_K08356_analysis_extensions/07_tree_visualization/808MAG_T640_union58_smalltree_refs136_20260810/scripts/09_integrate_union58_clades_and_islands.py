#!/usr/bin/env python3
"""Integrate tree/HMM evidence with +/-10 CDS neighborhoods for 58 candidates."""

import csv
from collections import defaultdict
from pathlib import Path


OUT = Path(__file__).resolve().parent.parent
PLACEMENT = OUT / "04_tables/candidate58_tree_placement_summary.tsv"
NEIGH = OUT / "04_neighborhood/union58_neighborhood_long.tsv"
HITS = OUT / "04_neighborhood/union58_joint_reference_hits.raw.tsv"
OLD_EVIDENCE = OUT / "04_tables/source_old49_strict_qc_evidence.tsv"
OLD_ORFS = OUT / "04_neighborhood/source_old48_gene_neighborhood_ORF_table.csv"


CLASS_MAP = {
    "canonical AioA-associated": "Canonical AioA clade",
    "unknown/uncertain DMSOR": "Unknown DMSOR clade",
    "IdrA-associated phylogenetic lineage": "IdrA-associated clade",
    "synteny-supported IdrA/DIRM-like neighborhood-associated lineage": "IdrA clade",
}
CLASS_COLORS = {
    "Canonical AioA clade": "#009E73",
    "Unknown DMSOR clade": "#D55E00",
    "IdrA-associated clade": "#0072B2",
    "IdrA clade": "#7B3294",
}


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def ref_family(subject_id):
    if subject_id.startswith(("IdrP_like", "IdrP1", "IdrP2")):
        return "P_like"
    if subject_id.startswith("AioB"):
        return "canonical_AioB"
    if subject_id.startswith("IdrB"):
        return "IdrB_related"
    return "other"


def accepted(fields):
    fam = ref_family(fields[1])
    evalue, bitscore = float(fields[10]), float(fields[11])
    length, qlen, slen = float(fields[3]), float(fields[12]), float(fields[13])
    qcov = 100 * length / qlen
    scov = 100 * length / slen
    if fam in {"IdrB_related", "canonical_AioB"}:
        ok = evalue <= 1e-5 and bitscore >= 50 and qcov >= 50 and scov >= 35
    elif fam == "P_like":
        ok = evalue <= 1e-10 and bitscore >= 80 and qcov >= 45 and scov >= 35
    else:
        ok = False
    return ok, fam, bitscore, qcov, scov, evalue


placement = read_tsv(PLACEMENT)
neighborhood = read_tsv(NEIGH)
old_evidence = {row["Sequence_ID"]: row for row in read_tsv(OLD_EVIDENCE)}

old_orf = {}
with OLD_ORFS.open(newline="") as handle:
    for row in csv.DictReader(handle):
        old_orf[(row["target_gene_id"], row["neighbor_orf_id"])] = row

best_hits = {}
with HITS.open(newline="") as handle:
    for fields in csv.reader(handle, delimiter="\t"):
        if len(fields) < 14:
            continue
        ok, family, bitscore, qcov, scov, evalue = accepted(fields)
        if not ok:
            continue
        key = (fields[0], family)
        if key not in best_hits or bitscore > best_hits[key]["bitscore"]:
            best_hits[key] = {
                "subject": fields[1], "family": family, "bitscore": bitscore,
                "qcov_pct": qcov, "scov_pct": scov, "evalue": evalue,
            }

rows_by_target = defaultdict(list)
for row in neighborhood:
    rows_by_target[row["target_id"]].append(row)

evidence_rows = []
for row in placement:
    seq_id = row["protein_ID"]
    neigh_rows = rows_by_target[seq_id]
    support = []
    for gene in neigh_rows:
        gene_id = gene["neighbor_raw_id"]
        for family in ("IdrB_related", "canonical_AioB", "P_like"):
            hit = best_hits.get((gene_id, family))
            if hit:
                support.append((gene, hit))
    b_hits = [(gene, hit) for gene, hit in support if hit["family"] in {"IdrB_related", "canonical_AioB"}]
    p_hits = [(gene, hit) for gene, hit in support if hit["family"] == "P_like"]
    b_hits.sort(key=lambda item: -item[1]["bitscore"])
    p_best = {}
    for gene, hit in p_hits:
        p_best.setdefault(gene["neighbor_raw_id"], (gene, hit))
    p_hits = sorted(p_best.values(), key=lambda item: int(item[0]["relative_position"]))

    old = old_evidence.get(seq_id)
    if old:
        final_class = CLASS_MAP[old["Final_class"]]
        confidence = old["Confidence"]
        rule = "retained from strict-QC 49-protein evidence table"
        strict = old["strict_synteny_QC_pass"]
        neighborhood_class = old["neighborhood_class"]
    else:
        strict_complete = bool(b_hits) and len(p_hits) >= 2
        if row["nearest_reference_group"] == "Unknown-clade reference":
            final_class = "Unknown DMSOR clade"
            confidence = "Medium"
            rule = "new T640 candidate; unknown-clade tree placement; neighborhood shown but not used as final functional proof"
        elif row["nearest_reference_group"] == "AioA reference":
            final_class = "Canonical AioA clade"
            confidence = "Medium"
            rule = "new T640 candidate; AioA-reference tree placement"
        elif strict_complete:
            final_class = "IdrA clade"
            confidence = "Medium"
            rule = "new T640 candidate; IdrA-reference placement plus B-related and >=2 P-like neighbors"
        else:
            final_class = "IdrA-associated clade"
            confidence = "Medium"
            rule = "new T640 candidate; IdrA-reference placement without complete B+2P neighborhood support"
        strict = "yes" if strict_complete else "no"
        if strict_complete:
            neighborhood_class = "complete DIRM-like neighborhood"
        elif b_hits and p_hits:
            neighborhood_class = "partial DIRM-like neighborhood"
        elif b_hits and b_hits[0][1]["family"] == "canonical_AioB":
            neighborhood_class = "Aio-like neighborhood"
        else:
            neighborhood_class = "unresolved/truncated neighborhood"

    evidence_rows.append({
        **row,
        "Final_class": final_class,
        "class_color": CLASS_COLORS[final_class],
        "Confidence": confidence,
        "classification_rule": rule,
        "strict_synteny_QC_pass": strict,
        "neighborhood_class": neighborhood_class,
        "B_related_gene": b_hits[0][0]["neighbor_raw_id"] if b_hits else "",
        "B_related_family": b_hits[0][1]["family"] if b_hits else "",
        "P_like_count": len(p_hits),
        "P_like_genes": ";".join(gene["neighbor_raw_id"] for gene, _ in p_hits),
    })

evidence_fields = list(placement[0]) + [
    "Final_class", "class_color", "Confidence", "classification_rule",
    "strict_synteny_QC_pass", "neighborhood_class", "B_related_gene",
    "B_related_family", "P_like_count", "P_like_genes",
]
write_tsv(OUT / "04_tables/candidate58_four_clade_integrated_evidence.tsv", evidence_rows, evidence_fields)

evidence_by_id = {row["protein_ID"]: row for row in evidence_rows}
manifest = []
for target_id, genes in rows_by_target.items():
    target = next(row for row in genes if int(row["relative_position"]) == 0)
    target_mid = (int(target["start"]) + int(target["end"])) / 2
    for gene in sorted(genes, key=lambda item: int(item["start"])):
        gene_id = gene["neighbor_raw_id"]
        old = old_orf.get((target_id, gene_id), {})
        if int(gene["relative_position"]) == 0:
            plot_label = "AioA/IdrA-related A"
        elif (gene_id, "IdrB_related") in best_hits or (gene_id, "canonical_AioB") in best_hits:
            plot_label = "B-related small subunit"
        elif (gene_id, "P_like") in best_hits:
            plot_label = "P-like"
        elif old.get("KO") or old.get("product"):
            plot_label = "other annotated CDS"
        else:
            plot_label = "function unknown"
        manifest.append({
            "tree_label": next(row["prefixed_target_ID"] for row in evidence_rows if row["protein_ID"] == target_id),
            "candidate_id": target_id,
            "MAG_ID": evidence_by_id[target_id]["MAG_ID"],
            "Final_class": evidence_by_id[target_id]["Final_class"],
            "gene_id": gene_id,
            "relative_position": gene["relative_position"],
            "relative_bp_start": f'{int(gene["start"]) - target_mid:.1f}',
            "relative_bp_end": f'{int(gene["end"]) - target_mid:.1f}',
            "strand": gene["strand"],
            "plot_label": plot_label,
            "KO": old.get("KO", ""),
            "annotation": old.get("product", ""),
        })

manifest_fields = list(manifest[0])
write_tsv(OUT / "04_neighborhood/union58_gene_island_manifest.tsv", manifest, manifest_fields)

counts = defaultdict(int)
for row in evidence_rows:
    counts[row["Final_class"]] += 1
print("class_counts=" + "; ".join(f"{key}:{counts[key]}" for key in CLASS_COLORS))
print(f"evidence_rows={len(evidence_rows)}")
print(f"manifest_rows={len(manifest)}")
