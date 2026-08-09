#!/usr/bin/env python3
import csv
import re
from pathlib import Path

import pandas as pd


BASE = Path("/Users/catherine/Downloads/MAG-bin-tpm")
OUTDIR = BASE / "derepMAG_abundance_stats"
TREE_DIR = Path("/Users/catherine/Downloads/aioA 蛋白序列建树/MAG 补充 idra-tree/latest_MAG_tree_beautify_files")
SOURCE_XLSX = TREE_DIR / "MAG_tree_clade_statistics_like_tree11.xlsx"
SOURCE_COLORSTRIP = TREE_DIR / "itol_COLORSTRIP_own_MAG_habitat.txt"

COLORS = {
    "IS": "#009E73",
    "AS": "#CC79A7",
    "ES": "#E69F00",
    "NS": "#0072B2",
}


def infer_habitat(sample_or_mag):
    sample = normalize_mag(sample_or_mag)
    sample = sample.split("_bin", 1)[0]
    if sample.startswith("SY"):
        return "IS"
    if sample.startswith("SQ_"):
        return "AS"
    if sample.startswith(("S13", "S14", "S15")):
        return "ES"
    if sample.startswith("S3"):
        return "NS"
    if sample.startswith(("S1", "S2", "S4")):
        return "AS"
    if sample.startswith(("C", "ES_")):
        return "ES"
    if sample.startswith(("NS_", "R2111_")):
        return "NS"
    return "Unknown"


def sample_group(sample_or_mag):
    sample = normalize_mag(sample_or_mag).split("_bin", 1)[0]
    if sample.startswith("R2111_"):
        return "_".join(sample.split("_")[:2])
    return sample.split("_", 1)[0] if sample.startswith("S") and not sample.startswith("SQ_") else sample.rsplit("-", 2)[0]


def normalize_mag(mag_id):
    s = str(mag_id).strip()
    s = re.sub(r"\s+", "_", s)
    if "-k141_" in s:
        s = s.split("-k141_", 1)[0]
    elif "-k141" in s:
        s = s.split("-k141", 1)[0]
    return s


def matrix_mags():
    with (BASE / "merged_tpm_matrix.txt").open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        return {row[0] for row in reader if row}


def write_branch_map():
    OUTDIR.mkdir(exist_ok=True)
    mags = matrix_mags()
    df = pd.read_excel(SOURCE_XLSX, sheet_name="MAG detail")
    rows = []
    missing = []
    changed = []
    for _, r in df.iterrows():
        clade = str(r["assigned_clade"]).strip()
        original = str(r["mag_id"]).strip()
        if clade == "nan" or original == "nan":
            continue
        mag = normalize_mag(original)
        corrected = infer_habitat(mag)
        source_habitat = str(r.get("habitat", "")).strip()
        group = sample_group(mag)
        rows.append([mag, clade, original, source_habitat, corrected, group])
        if mag not in mags:
            missing.append([mag, clade, original])
        if source_habitat and source_habitat != corrected:
            changed.append([original, source_habitat, corrected])

    with (OUTDIR / "k08356_branch_feature_map.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "MAG",
            "Feature",
            "Original_MAG_ID",
            "Source_habitat_in_tree_table",
            "Corrected_habitat",
            "Corrected_sample_group",
        ])
        w.writerows(rows)

    with (OUTDIR / "tree_habitat_corrections.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Original_MAG_ID", "Tree_file_habitat", "Corrected_habitat"])
        w.writerows(changed)

    if missing:
        with (OUTDIR / "k08356_branch_map_missing_mags.tsv").open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["MAG", "Feature", "Original_MAG_ID"])
            w.writerows(missing)


def write_corrected_colorstrip():
    lines = SOURCE_COLORSTRIP.read_text().splitlines()
    out_lines = []
    in_data = False
    audit_rows = []
    for line in lines:
        if line == "DATA":
            in_data = True
            out_lines.append(line)
            continue
        if not in_data or not line.strip():
            out_lines.append(line)
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            label, old_color, old_habitat = parts[:3]
            corrected = infer_habitat(label)
            corrected_color = COLORS.get(corrected, old_color)
            out_lines.append("\t".join([label, corrected_color, corrected] + parts[3:]))
            audit_rows.append([label, old_habitat, corrected, "changed" if old_habitat != corrected else "ok"])
        else:
            out_lines.append(line)

    (OUTDIR / "itol_COLORSTRIP_own_MAG_habitat.corrected.txt").write_text("\n".join(out_lines) + "\n")
    with (OUTDIR / "itol_habitat_colorstrip_audit.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Tree_label", "Original_habitat", "Corrected_habitat", "Status"])
        w.writerows(audit_rows)


def main():
    write_branch_map()
    write_corrected_colorstrip()
    print("Prepared corrected branch map and tree habitat audit in", OUTDIR)


if __name__ == "__main__":
    main()
