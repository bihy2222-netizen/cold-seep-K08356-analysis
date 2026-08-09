#!/usr/bin/env python3
import csv
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


BASE = Path("/Users/catherine/Downloads/MAG-bin-tpm")
SOURCE = Path("/Users/catherine/Downloads/aioA 蛋白序列建树/MAG 补充 idra-tree/latest_MAG_tree_beautify_files")
DEST = BASE / "checked_corrected_tree_beautify_files_20260708"
ORIG = DEST / "ORIGINAL_BACKUP"
CORR = DEST / "CORRECTED_USE_THESE"
AUDIT = DEST / "AUDIT_TABLES"

COLORS = {
    "IS": "#009E73",
    "AS": "#CC79A7",
    "ES": "#E69F00",
    "NS": "#0072B2",
}

FORMAL_FILES = [
    "itol_TREE_COLORS_highlight_own_MAG_and_refs.txt",
    "latest_MAG_tree_for_iTOL.treefile",
    "itol_COLORSTRIP_own_MAG_habitat.txt",
    "MAG_tree_clade_lists_preview.png",
    "annotation_summary.tsv",
    "build_MAG_tree_clade_workbook.mjs",
    "classify_mag_tree_clades.py",
    "MAG_tree_clade_statistics_like_tree11.xlsx",
    "MAG_tree_clade_detail.csv",
    "README_iTOL_upload_order.txt",
    "itol_COLORSTRIP_sequence_type.txt",
    "make_latest_mag_tree_itol_annotations.py",
    "itol_COLORSTRIP_reference_sources.txt",
    "MAG_tree_clade_lists_like_tree11.tsv",
    "MAG_tree_clade_statistics.json",
]


def normalize_mag(label):
    s = str(label).strip()
    s = re.sub(r"\s+", "_", s)
    if "-k141_" in s:
        return s.split("-k141_", 1)[0]
    if "-k141" in s:
        return s.split("-k141", 1)[0]
    return s.split("\t", 1)[0]


def infer_habitat(label):
    mag = normalize_mag(label)
    sample = mag.split("_bin", 1)[0]
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


def corrected_group(label):
    mag = normalize_mag(label)
    sample = mag.split("_bin", 1)[0]
    if sample.startswith("R2111_"):
        return "_".join(sample.split("_")[:2])
    if sample.startswith("SQ_"):
        return "_".join(sample.split("_")[:2])
    if sample.startswith("SY"):
        return re.match(r"^(SY\d+[A-Z]+)", sample).group(1) if re.match(r"^(SY\d+[A-Z]+)", sample) else sample
    if sample.startswith(("S", "C")):
        return sample.split("_", 1)[0]
    if sample.startswith("ES_"):
        return "_".join(sample.split("_")[:2])
    if sample.startswith("NS_"):
        return "NS"
    return sample


def display(original_id):
    return f"{original_id} | {infer_habitat(original_id)} | {corrected_group(original_id)}"


def copy_originals():
    ORIG.mkdir(parents=True, exist_ok=True)
    for name in FORMAL_FILES:
        src = SOURCE / name
        if src.exists():
            shutil.copy2(src, ORIG / name)


def copy_passthrough_corrected():
    CORR.mkdir(parents=True, exist_ok=True)
    passthrough = [
        "latest_MAG_tree_for_iTOL.treefile",
        "itol_COLORSTRIP_sequence_type.txt",
        "itol_COLORSTRIP_reference_sources.txt",
        "MAG_tree_clade_lists_preview.png",
        "build_MAG_tree_clade_workbook.mjs",
        "classify_mag_tree_clades.py",
        "make_latest_mag_tree_itol_annotations.py",
    ]
    for name in passthrough:
        src = SOURCE / name
        if src.exists():
            shutil.copy2(src, CORR / name)


def correct_colorstrip():
    src = SOURCE / "itol_COLORSTRIP_own_MAG_habitat.txt"
    out = CORR / "itol_COLORSTRIP_own_MAG_habitat.txt"
    audit = []
    lines = src.read_text().splitlines()
    new = []
    in_data = False
    for line in lines:
        if line == "DATA":
            in_data = True
            new.append(line)
            continue
        if not in_data or not line.strip():
            new.append(line)
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            label, old_color, old_habitat = parts[:3]
            habitat = infer_habitat(label)
            color = COLORS.get(habitat, old_color)
            new.append("\t".join([label, color, habitat] + parts[3:]))
            audit.append([label, old_habitat, habitat, old_color, color, "changed" if old_habitat != habitat or old_color != color else "ok"])
        else:
            new.append(line)
    out.write_text("\n".join(new) + "\n")
    return audit


def correct_tree_colors():
    src = SOURCE / "itol_TREE_COLORS_highlight_own_MAG_and_refs.txt"
    out = CORR / "itol_TREE_COLORS_highlight_own_MAG_and_refs.txt"
    audit = []
    new = []
    for line in src.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 4 and "_bin" in parts[0] and "-k141" in parts[0]:
            label = parts[0]
            habitat = infer_habitat(label)
            color = COLORS.get(habitat, parts[2])
            old_color = parts[2]
            parts[2] = color
            audit.append([label, parts[1], old_color, color, habitat, "changed" if old_color != color else "ok"])
            new.append("\t".join(parts))
        else:
            new.append(line)
    out.write_text("\n".join(new) + "\n")
    return audit


def correct_detail_csv():
    df = pd.read_csv(SOURCE / "MAG_tree_clade_detail.csv")
    df["original_habitat"] = df["habitat"]
    df["habitat"] = df["mag_id"].map(infer_habitat)
    df["sample_group"] = df["mag_id"].map(corrected_group)
    df["display"] = df["mag_id"].map(display)
    df.to_csv(CORR / "MAG_tree_clade_detail.csv", index=False)
    changed = df.loc[df["original_habitat"] != df["habitat"], ["mag_id", "original_habitat", "habitat", "assigned_clade"]]
    changed.to_csv(AUDIT / "tree_habitat_corrections.tsv", sep="\t", index=False)
    return df, changed


def write_clade_lists(df):
    order = [
        "IdrA-associated",
        "canonical aioA-associated",
        "aioA-like-associated",
        "unknown AioA-like / uncertain DMSOR",
        "other DMSOR-family reference-associated",
    ]
    headers = []
    cols = []
    for clade in order:
        values = df.loc[df["assigned_clade"] == clade, "display"].tolist()
        headers.append(f"{clade} (n={len(values)})")
        cols.append(values)
    max_len = max([len(c) for c in cols] + [0])
    rows = []
    for i in range(max_len):
        rows.append([col[i] if i < len(col) else "" for col in cols])
    with (CORR / "MAG_tree_clade_lists_like_tree11.tsv").open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(rows)
    return headers, rows


def write_annotation_summary(df):
    counts = Counter(df["habitat"])
    ref_counts = {}
    source_summary = SOURCE / "annotation_summary.tsv"
    if source_summary.exists():
        with source_summary.open(newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row["category"].startswith("reference:"):
                    ref_counts[row["category"]] = int(row["count"])
    with (CORR / "annotation_summary.tsv").open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["category", "count"])
        for habitat in ["AS", "ES", "IS", "NS"]:
            writer.writerow([f"habitat:{habitat}", counts.get(habitat, 0)])
        for key in sorted(ref_counts):
            writer.writerow([key, ref_counts[key]])


def write_workbook(df, headers, clade_rows):
    summary = df["assigned_clade"].value_counts().rename_axis("assigned_clade").reset_index(name="n")
    order = {
        "IdrA-associated": 0,
        "canonical aioA-associated": 1,
        "aioA-like-associated": 2,
        "unknown AioA-like / uncertain DMSOR": 3,
        "other DMSOR-family reference-associated": 4,
    }
    summary["_order"] = summary["assigned_clade"].map(order).fillna(99)
    summary = summary.sort_values("_order").drop(columns="_order")
    summary = pd.concat([summary, pd.DataFrame([{"assigned_clade": "Total MAG/sample sequences", "n": int(summary["n"].sum())}])], ignore_index=True)
    habitat_summary = df.groupby(["assigned_clade", "habitat"], as_index=False).size().rename(columns={"size": "n"})
    habitat_summary["_order"] = habitat_summary["assigned_clade"].map(order).fillna(99)
    habitat_summary = habitat_summary.sort_values(["_order", "habitat"]).drop(columns="_order")
    readme = pd.DataFrame({
        "Field": [
            "assigned_clade",
            "nearest_reference",
            "nearest_3_references",
            "nearest_5_clade_votes",
            "Correction note",
        ],
        "Description": [
            "MAG/sample sequence clade assigned by the nearest non-MAG reference in the IQ-TREE branch-length tree.",
            "Closest published/reference sequence by tree distance.",
            "Top 3 nearest references for audit; tree tip names are unchanged.",
            "Reference-clade composition among the 5 nearest references.",
            "Habitat labels corrected by study design: S13 and S15 are ES; S3/R2111 are NS; SQ/S1/S2/S4 are AS; SY are IS.",
        ],
    })
    clade_df = pd.DataFrame(clade_rows, columns=headers)
    detail = df.drop(columns=["original_habitat"], errors="ignore")
    with pd.ExcelWriter(CORR / "MAG_tree_clade_statistics_like_tree11.xlsx", engine="openpyxl") as writer:
        clade_df.to_excel(writer, sheet_name="Clade lists", index=False)
        detail.to_excel(writer, sheet_name="MAG detail", index=False)
        summary.to_excel(writer, sheet_name="Summary", index=False)
        habitat_summary.to_excel(writer, sheet_name="Habitat summary", index=False)
        readme.to_excel(writer, sheet_name="README", index=False)


def write_corrected_json(df, headers, clade_rows):
    counts = Counter(df["assigned_clade"])
    habitat_counts = (
        df.groupby(["assigned_clade", "habitat"], as_index=False)
        .size()
        .rename(columns={"size": "n"})
        .to_dict(orient="records")
    )
    detail = df.drop(columns=["original_habitat"], errors="ignore").to_dict(orient="records")
    data = {
        "matrix": [headers] + clade_rows,
        "detail": detail,
        "counts": dict(counts),
        "habitat_counts": habitat_counts,
        "clade_order": [
            "IdrA-associated",
            "canonical aioA-associated",
            "aioA-like-associated",
            "unknown AioA-like / uncertain DMSOR",
            "other DMSOR-family reference-associated",
        ],
        "correction_note": "Habitat labels corrected by study design on 2026-07-08; S13/S15 changed from AS to ES.",
        "corrected_habitat_counts": dict(Counter(df["habitat"])),
    }
    (CORR / "MAG_tree_clade_statistics.json").write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def write_readme(color_audit, tree_color_audit, changed):
    changed_lines = ["Corrected habitat mismatches:"]
    if len(changed) == 0:
        changed_lines.append("- None")
    else:
        for _, row in changed.iterrows():
            changed_lines.append(f"- {row['mag_id']}: {row['original_habitat']} -> {row['habitat']}")
    text = "\n".join([
        "Checked/corrected MAG tree beautify files",
        "",
        "Use files in CORRECTED_USE_THESE for iTOL upload and downstream tables.",
        "ORIGINAL_BACKUP contains the original files copied from latest_MAG_tree_beautify_files.",
        "",
        "iTOL upload order:",
        "1. latest_MAG_tree_for_iTOL.treefile",
        "2. itol_COLORSTRIP_sequence_type.txt",
        "3. itol_COLORSTRIP_own_MAG_habitat.txt",
        "4. itol_COLORSTRIP_reference_sources.txt",
        "5. itol_TREE_COLORS_highlight_own_MAG_and_refs.txt",
        "",
        *changed_lines,
        "",
        "Audit tables:",
        "- AUDIT_TABLES/tree_habitat_corrections.tsv",
        "- AUDIT_TABLES/itol_habitat_colorstrip_audit.tsv",
        "- AUDIT_TABLES/itol_tree_colors_audit.tsv",
        "",
        f"Own MAG habitat colorstrip rows checked: {len(color_audit)}",
        f"Own MAG TREE_COLORS rows checked: {len(tree_color_audit)}",
    ])
    (DEST / "README_CHECKED_CORRECTED.txt").write_text(text + "\n")
    (CORR / "README_iTOL_upload_order.txt").write_text(text + "\n")


def write_audits(color_audit, tree_color_audit):
    AUDIT.mkdir(parents=True, exist_ok=True)
    with (AUDIT / "itol_habitat_colorstrip_audit.tsv").open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Tree_label", "Original_habitat", "Corrected_habitat", "Original_color", "Corrected_color", "Status"])
        writer.writerows(color_audit)
    with (AUDIT / "itol_tree_colors_audit.tsv").open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Tree_label", "Element", "Original_color", "Corrected_color", "Corrected_habitat", "Status"])
        writer.writerows(tree_color_audit)


def copy_report_files():
    report_files = [
        "K08356_derepMAG_significance_report.xlsx",
        "sample_habitat_map.tsv",
        "sample_group_summary.tsv",
        "k08356_branch_feature_map.tsv",
    ]
    report_dir = DEST / "ABUNDANCE_REPORT"
    report_dir.mkdir(parents=True, exist_ok=True)
    for name in report_files:
        src = BASE / "derepMAG_abundance_stats" / name
        if src.exists():
            shutil.copy2(src, report_dir / name)


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    ORIG.mkdir(parents=True, exist_ok=True)
    CORR.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    copy_originals()
    copy_passthrough_corrected()
    color_audit = correct_colorstrip()
    tree_color_audit = correct_tree_colors()
    df, changed = correct_detail_csv()
    headers, clade_rows = write_clade_lists(df)
    write_annotation_summary(df)
    write_workbook(df, headers, clade_rows)
    write_corrected_json(df, headers, clade_rows)
    write_audits(color_audit, tree_color_audit)
    write_readme(color_audit, tree_color_audit, changed)
    copy_report_files()
    print(DEST)


if __name__ == "__main__":
    main()
