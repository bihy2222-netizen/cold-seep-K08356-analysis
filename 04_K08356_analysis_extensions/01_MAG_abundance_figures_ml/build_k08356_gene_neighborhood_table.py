#!/usr/bin/env python3
import csv
import os
import re
from collections import defaultdict

BASE = "/Users/catherine/Downloads/MAG-bin-tpm/K08356_MAG_gene_island_arrows_20260708"
DATA_DIR = os.path.join(BASE, "sample_cds_gff")
BRANCH_MAP = "/Users/catherine/Downloads/MAG-bin-tpm/derepMAG_abundance_stats/k08356_branch_feature_map.tsv"
OUT_COORDS = os.path.join(BASE, "K08356_four_branch_gene_neighborhood_coordinates.csv")
OUT_STATUS = os.path.join(BASE, "K08356_four_branch_gene_neighborhood_status.csv")
FLANK = 10


def sample_from_original(original_id):
    return re.sub(r"_bin\d+-k141_.*$", "", original_id)


def target_short_id(original_id):
    m = re.search(r"(k141_\d+_\d+)$", original_id)
    return m.group(1) if m else None


def contig_from_short(short_id):
    return short_id.rsplit("_", 1)[0]


def parse_faa_headers(path, contig):
    records = []
    if not os.path.exists(path):
        return records
    with open(path, "r", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            # >k141_14765_13 # 14352 # 16832 # 1 # ID=...
            m = re.match(r"^>(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#", line)
            if not m:
                continue
            prot, start, end, strand = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
            if contig_from_short(prot) != contig:
                continue
            records.append({
                "Protein_short": prot,
                "Contig": contig,
                "ORF_number": int(prot.rsplit("_", 1)[1]),
                "Start": start,
                "End": end,
                "Direction": strand,
            })
    records.sort(key=lambda r: (r["Start"], r["End"], r["ORF_number"]))
    return records


def classify_annotation(ko, desc, is_target, feature):
    text = f"{ko} {desc}".lower()
    if is_target:
        if feature == "IdrA-associated":
            return "Target IdrA-like K08356", "IdrA"
        if feature == "canonical aioA-associated":
            return "Target canonical aioA", "AioA"
        if feature == "aioA-like-associated":
            return "Target aioA-like K08356", "AioA-like"
        return "Target uncertain AioA-like", "AioA-like?"
    if ko == "K08355":
        return "Rieske/AioB", "AioB"
    if "rieske" in text or "2fe-2s" in text or "iron-sulfur" in text:
        return "Rieske/AioB", "Rieske"
    if "cytochrome" in text or "peroxidase" in text:
        return "Cytochrome/peroxidase", "cyt"
    if "arsen" in text or "acr3" in text or "arsb" in text or "arsc" in text or "arsh" in text:
        return "Arsenic resistance", "ars"
    if "transporter" in text or "permease" in text or "efflux" in text or "pump" in text or "pst" in text or "mfs" in text:
        return "Transporter", "transporter"
    if "sox" in text or "sulfur oxidation" in text or "sulphur oxidation" in text or "thiosulfate" in text:
        return "Sox/sulfur oxidation", "sox"
    if "nir" in text or "nitrite" in text or "nitrate" in text or "denitrification" in text:
        return "Nitrogen metabolism", "nir"
    if "ferredoxin" in text or "oxidoreductase" in text or "dehydrogenase" in text:
        return "Redox/oxidoreductase", "Redox"
    if desc:
        return "Annotated other", ""
    return "Hypothetical/unknown", ""


def short_desc(desc):
    desc = desc.strip().strip('"')
    desc = re.sub(r"\s+\[EC:.*?\]", "", desc)
    desc = re.sub(r"\s+", " ", desc)
    return desc


targets = []
with open(BRANCH_MAP, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        original = row["Original_MAG_ID"]
        sample = sample_from_original(original)
        short_id = target_short_id(original)
        targets.append({
            "MAG": row["MAG"],
            "Feature": row["Feature"],
            "Original_MAG_ID": original,
            "Corrected_habitat": row["Corrected_habitat"],
            "Corrected_sample_group": row["Corrected_sample_group"],
            "Sample": sample,
            "Target_short": short_id,
            "Target_contig": contig_from_short(short_id) if short_id else None,
        })

rows = []
status = []
ids_needed_by_sample = defaultdict(set)
pending_neighborhoods = []

for t in targets:
    sample = t["Sample"]
    faa = os.path.join(DATA_DIR, f"{sample}.cdhit.kofamscan.best.contig.cds.faa")
    records = parse_faa_headers(faa, t["Target_contig"])
    target_indices = [i for i, r in enumerate(records) if r["Protein_short"] == t["Target_short"]]
    if not records or not target_indices:
        status.append({**t, "Status": "target_not_found", "n_contig_orfs": len(records), "n_neighborhood_orfs": 0})
        continue
    idx = target_indices[0]
    lo = max(0, idx - FLANK)
    hi = min(len(records), idx + FLANK + 1)
    subset = records[lo:hi]
    status.append({**t, "Status": "ok", "n_contig_orfs": len(records), "n_neighborhood_orfs": len(subset)})
    for local_order, r in enumerate(subset, start=1):
        ids_needed_by_sample[sample].add(r["Protein_short"])
        pending_neighborhoods.append((t, r, local_order, idx - lo + 1))

annotations = {}
for sample, needed in ids_needed_by_sample.items():
    best = os.path.join(DATA_DIR, f"{sample}.cdhit.kofamscan.best.txt")
    if not os.path.exists(best):
        continue
    with open(best, "r", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            prot = parts[1]
            if prot in needed:
                annotations[(sample, prot)] = {
                    "KO": parts[2],
                    "Score": parts[3],
                    "E_value": parts[5],
                    "Description": parts[6].strip().strip('"'),
                }

for t, r, local_order, target_local_order in pending_neighborhoods:
    ann = annotations.get((t["Sample"], r["Protein_short"]), {})
    is_target = r["Protein_short"] == t["Target_short"]
    gene_class, gene_label = classify_annotation(ann.get("KO", ""), ann.get("Description", ""), is_target, t["Feature"])
    if not gene_label and is_target:
        gene_label = "K08356"
    molecule = f"{t['Original_MAG_ID']} | {t['Feature']} | {t['Corrected_habitat']}"
    rows.append({
        **t,
        **r,
        "Molecule": molecule,
        "Local_order": local_order,
        "Target_local_order": target_local_order,
        "Is_target": "yes" if is_target else "no",
        "KO": ann.get("KO", ""),
        "Description": short_desc(ann.get("Description", "")),
        "Gene_class": gene_class,
        "Gene_label": gene_label,
        "Mid": (r["Start"] + r["End"]) / 2,
    })

fieldnames = [
    "MAG", "Feature", "Original_MAG_ID", "Corrected_habitat", "Corrected_sample_group", "Sample",
    "Target_short", "Target_contig", "Molecule", "Protein_short", "Contig", "ORF_number",
    "Start", "End", "Mid", "Direction", "Local_order", "Target_local_order", "Is_target",
    "KO", "Description", "Gene_class", "Gene_label"
]
with open(OUT_COORDS, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

status_fields = [
    "MAG", "Feature", "Original_MAG_ID", "Corrected_habitat", "Corrected_sample_group",
    "Sample", "Target_short", "Target_contig", "Status", "n_contig_orfs", "n_neighborhood_orfs"
]
with open(OUT_STATUS, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=status_fields)
    writer.writeheader()
    writer.writerows(status)

print(f"Wrote {len(rows)} neighborhood ORFs")
print(f"Wrote status for {len(status)} targets")
