#!/usr/bin/env python3
"""Overlay the four evidence-defined MAG clades on the small-tree SVG."""

import argparse
import csv
import html
import re
from collections import Counter
from pathlib import Path


CLASS_ORDER = [
    "canonical AioA-associated",
    "synteny-supported IdrA/DIRM-like neighborhood-associated lineage",
    "IdrA-associated phylogenetic lineage",
    "unknown/uncertain DMSOR",
]

CLASS_COLORS = {
    "canonical AioA-associated": "#009E73",
    "synteny-supported IdrA/DIRM-like neighborhood-associated lineage": "#7B3294",
    "IdrA-associated phylogenetic lineage": "#0072B2",
    "unknown/uncertain DMSOR": "#D55E00",
}

CLASS_LABELS = {
    "canonical AioA-associated": "Canonical AioA",
    "synteny-supported IdrA/DIRM-like neighborhood-associated lineage": "DIRM-synteny IdrA",
    "IdrA-associated phylogenetic lineage": "IdrA phylogenetic",
    "unknown/uncertain DMSOR": "Uncertain DMSOR",
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


def esc(value):
    return html.escape(str(value), quote=True)


def candidate_aliases(tree_label):
    candidate = tree_label.split("|", 1)[1] if "|" in tree_label else tree_label
    return {tree_label, candidate, candidate.replace("_bin", " bin", 1)}


def add_inline_attribute(attributes, name, value):
    pattern = re.compile(rf'\s{name}="[^"]*"')
    replacement = f' {name}="{esc(value)}"'
    if pattern.search(attributes):
        return pattern.sub(replacement, attributes, count=1)
    return attributes + replacement


def make_legend(counts):
    x0 = 350
    y_title = 20
    y = 38
    chunks = [
        '<g id="codex-four-clade-legend" font-family="Arial, Helvetica, sans-serif">',
        f'<text x="{x0}" y="{y_title}" font-size="11" font-weight="700" fill="#111827">Evidence-defined MAG clades</text>',
    ]
    x = x0
    widths = [160, 190, 180, 165]
    for class_name, width in zip(CLASS_ORDER, widths):
        color = CLASS_COLORS[class_name]
        label = f"{CLASS_LABELS[class_name]} (n={counts[class_name]})"
        chunks.append(
            f'<rect x="{x}" y="{y - 9}" width="13" height="10" fill="{color}" '
            'stroke="#FFFFFF" stroke-width="0.6"/>'
        )
        chunks.append(
            f'<text x="{x + 18}" y="{y}" font-size="9" font-weight="400" fill="#374151">{esc(label)}</text>'
        )
        x += width
    chunks.append("</g>")
    return "\n".join(chunks)


def overlay_classes(svg, evidence_rows):
    by_tree_label = {}
    aliases = {}
    for row in evidence_rows:
        tree_label = f'{row["Habitat"]}|{row["Sequence_ID"]}'
        row = dict(row)
        row["tree_label"] = tree_label
        by_tree_label[tree_label] = row
        for alias in candidate_aliases(tree_label):
            aliases[alias] = tree_label

    pattern = re.compile(
        r'<text\b(?P<attrs>[^>]*transform="translate\((?P<x>[-0-9.]+)[ ,]+(?P<y>[-0-9.]+)\)"[^>]*)>'
        r'(?P<body><tspan[^>]*>(?P<label>.*?)</tspan>)</text>',
        re.S,
    )
    matched = {}

    def replace_text(match):
        visible = re.sub(r"<.*?>", "", match.group("label")).strip()
        visible_candidate = visible.split("|", 1)[1] if "|" in visible else visible
        tree_label = aliases.get(visible) or aliases.get(visible_candidate)
        if not tree_label:
            return match.group(0)
        row = by_tree_label[tree_label]
        class_name = row["Final_class"]
        color = CLASS_COLORS[class_name]
        attributes = match.group("attrs")
        attributes = add_inline_attribute(attributes, "fill", color)
        attributes = add_inline_attribute(attributes, "font-weight", "700")
        attributes = add_inline_attribute(attributes, "data-final-class", class_name)
        matched[tree_label] = {
            "x": float(match.group("x")),
            "y": float(match.group("y")),
            "visible_label": visible,
            "source_svg_tree_label": visible,
        }
        normalized_body = match.group("body").replace(match.group("label"), esc(tree_label), 1)
        return f'<text{attributes}>{normalized_body}</text>'

    svg = pattern.sub(replace_text, svg)
    missing = sorted(set(by_tree_label) - set(matched))
    if missing:
        raise SystemExit("Evidence rows without SVG labels: " + ", ".join(missing))

    strip_x = 1169.0
    strips = ['<g id="codex-four-clade-strip">']
    for tree_label, pos in sorted(matched.items(), key=lambda item: item[1]["y"]):
        row = by_tree_label[tree_label]
        color = CLASS_COLORS[row["Final_class"]]
        y = pos["y"] - 8.2
        strips.append(
            f'<rect x="{strip_x:.1f}" y="{y:.2f}" width="13" height="12" rx="1" '
            f'fill="{color}" stroke="#FFFFFF" stroke-width="0.7" '
            f'data-sequence-id="{esc(row["Sequence_ID"])}" data-confidence="{esc(row["Confidence"])}"/>'
        )
    strips.append("</g>")

    counts = Counter(row["Final_class"] for row in evidence_rows)
    overlay = "\n".join(strips) + "\n" + make_legend(counts)
    svg = svg.replace("</svg>", overlay + "\n</svg>")
    return svg, by_tree_label, matched


def classification_rows(evidence_rows, matched):
    selected = [
        "Sequence_ID",
        "MAG_ID",
        "Habitat",
        "Sample_group",
        "Final_class",
        "Confidence",
        "Tree_clade",
        "nearest_reference",
        "nearest_reference_family",
        "Combined_score",
        "IdrA_score",
        "AioA_score",
        "IdrA_minus_AioA",
        "neighborhood_class",
        "strict_synteny_QC_pass",
        "strict_synteny_QC_fail_reasons",
        "gene_order",
    ]
    rows = []
    for row in evidence_rows:
        out = {field: row.get(field, "") for field in selected}
        out["tree_label"] = f'{row["Habitat"]}|{row["Sequence_ID"]}'
        out["source_svg_tree_label"] = matched[out["tree_label"]]["source_svg_tree_label"]
        if "|" not in out["source_svg_tree_label"]:
            out["habitat_label_status"] = "prefix_missing_in_source_svg"
        else:
            svg_habitat = out["source_svg_tree_label"].split("|", 1)[0]
            out["habitat_label_status"] = "match" if svg_habitat == row["Habitat"] else "corrected_to_evidence_table"
        out["clade_display"] = CLASS_LABELS[row["Final_class"]]
        out["clade_color"] = CLASS_COLORS[row["Final_class"]]
        rows.append(out)
    rows.sort(key=lambda row: (CLASS_ORDER.index(row["Final_class"]), row["Habitat"], row["Sequence_ID"]))
    fields = [
        "tree_label",
        "source_svg_tree_label",
        "habitat_label_status",
        "clade_display",
        "clade_color",
    ] + selected
    return rows, fields


def summary_rows(evidence_rows):
    rows = []
    for class_name in CLASS_ORDER:
        subset = [row for row in evidence_rows if row["Final_class"] == class_name]
        habitat_counts = Counter(row["Habitat"] for row in subset)
        confidence_counts = Counter(row["Confidence"] for row in subset)
        rows.append(
            {
                "Final_class": class_name,
                "clade_display": CLASS_LABELS[class_name],
                "clade_color": CLASS_COLORS[class_name],
                "n": len(subset),
                "IS": habitat_counts["IS"],
                "AS": habitat_counts["AS"],
                "ES": habitat_counts["ES"],
                "NS": habitat_counts["NS"],
                "High_confidence": confidence_counts["High"],
                "Medium_confidence": confidence_counts["Medium"],
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-svg", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output-svg", type=Path, required=True)
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    args = parser.parse_args()

    evidence_rows = read_tsv(args.evidence)
    unknown_classes = sorted(set(row["Final_class"] for row in evidence_rows) - set(CLASS_ORDER))
    if unknown_classes:
        raise SystemExit("Unexpected Final_class values: " + ", ".join(unknown_classes))

    svg = args.input_svg.read_text()
    svg, _, matched = overlay_classes(svg, evidence_rows)
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(svg)

    detail_rows, detail_fields = classification_rows(evidence_rows, matched)
    write_tsv(args.output_table, detail_rows, detail_fields)
    summary = summary_rows(evidence_rows)
    write_tsv(
        args.output_summary,
        summary,
        [
            "Final_class",
            "clade_display",
            "clade_color",
            "n",
            "IS",
            "AS",
            "ES",
            "NS",
            "High_confidence",
            "Medium_confidence",
        ],
    )

    manual_aliases = [
        f'{pos["visible_label"]} -> {tree_label}'
        for tree_label, pos in matched.items()
        if pos["visible_label"] != tree_label
    ]
    print(f"classified_svg_rows={len(matched)}")
    print("class_counts=" + "; ".join(
        f"{CLASS_LABELS[class_name]}:{sum(row['Final_class'] == class_name for row in evidence_rows)}"
        for class_name in CLASS_ORDER
    ))
    print("manual_label_aliases=" + ", ".join(manual_aliases))
    print(f"output_svg={args.output_svg}")


if __name__ == "__main__":
    main()
