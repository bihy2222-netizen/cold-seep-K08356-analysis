#!/usr/bin/env python3
"""Add aligned gene-island arrows to the right of the supplied Illustrator SVG tree."""

import argparse
import csv
import html
import math
import re
from collections import defaultdict
from pathlib import Path


GENE_COLORS = {
    "AioA/IdrA-related A": "#2166AC",
    "CD-search only A": "#2166AC",
    "IdrB/AioB-related small subunit": "#F58518",
    "canonical AioB": "#B279A2",
    "P-like": "#2CA25F",
    "function_unknown": "url(#geneIslandUnknownHatch)",
    "other annotated CDS": "#D7DCE2",
}

SHORT_LABELS = {
    "AioA/IdrA-related A": "A",
    "CD-search only A": "A",
    "IdrB/AioB-related small subunit": "B",
    "canonical AioB": "B",
    "P-like": "P",
}


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def esc(value):
    return html.escape(str(value), quote=True)


def candidate_aliases(tree_label):
    candidate = tree_label.split("|", 1)[1] if "|" in tree_label else tree_label
    aliases = {tree_label, candidate}
    aliases.add(candidate.replace("_bin", " bin", 1))
    return aliases


def parse_tip_positions(svg_text, canonical_labels):
    aliases = {}
    for label in canonical_labels:
        for alias in candidate_aliases(label):
            aliases[alias] = label

    pattern = re.compile(
        r'<text\b[^>]*transform="translate\(([-0-9.]+)[ ,]+([-0-9.]+)\)"[^>]*>'
        r'<tspan[^>]*>(.*?)</tspan></text>',
        re.S,
    )
    positions = {}
    source_labels = {}
    for match in pattern.finditer(svg_text):
        visible = re.sub(r"<.*?>", "", match.group(3)).strip()
        canonical = aliases.get(visible)
        if canonical:
            positions[canonical] = (float(match.group(1)), float(match.group(2)))
            source_labels[canonical] = visible
    return positions, source_labels


def augment_r2111(rows_by_label):
    label = "NS|R2111_N500_0-10_bin13-k141_4404301_3"
    rows = rows_by_label.get(label)
    if not rows:
        return
    if any(row.get("plot_label") in {"canonical AioB", "IdrB/AioB-related small subunit"} for row in rows):
        return
    rows.insert(
        0,
        {
            "tree_label": label,
            "candidate_id": "R2111_N500_0-10_bin13-k141_4404301_3",
            "gene_id": "R2111_N500_0-10_bin13-k141_4404301_2",
            "relative_bp_start": "-1827",
            "relative_bp_end": "-1237",
            "strand": "+",
            "plot_label": "canonical AioB",
            "annotation": "K08355 arsenite oxidase small subunit",
        },
    )


def arrow_points(x1, x2, y, height, strand):
    if x2 < x1:
        x1, x2 = x2, x1
    width = max(x2 - x1, 3.2)
    head = min(7.0, width * 0.42)
    top = y - height / 2
    bottom = y + height / 2
    if strand == "-":
        points = [(x1, y), (x1 + head, top), (x2, top), (x2, bottom), (x1 + head, bottom)]
    else:
        points = [(x1, top), (x2 - head, top), (x2, y), (x2 - head, bottom), (x1, bottom)]
    return " ".join(f"{x:.2f},{py:.2f}" for x, py in points)


def text(x, y, value, size=9, color="#111827", weight="400", anchor="start"):
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" '
        f'fill="{color}">{esc(value)}</text>'
    )


def build_overlay(rows_by_label, positions, panel_x, panel_width):
    displayed = sorted(positions, key=lambda label: positions[label][1])
    values = []
    for label in displayed:
        for row in rows_by_label[label]:
            try:
                values.extend((float(row["relative_bp_start"]), float(row["relative_bp_end"])))
            except (KeyError, TypeError, ValueError):
                continue

    extent = max([12500.0] + [abs(value) for value in values])
    extent = math.ceil(extent / 1000.0) * 1000.0
    center_x = panel_x + panel_width / 2

    def sx(relative_bp):
        return center_x + float(relative_bp) / (2 * extent) * panel_width

    y_min = min(positions[label][1] for label in displayed)
    y_max = max(positions[label][1] for label in displayed)
    parts = [
        '<g id="codex-gene-island-overlay">',
        '<defs>',
        '<pattern id="geneIslandUnknownHatch" width="3" height="3" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">',
        '<rect width="3" height="3" fill="#6B7280"/>',
        '<line x1="0" y1="0" x2="0" y2="3" stroke="#F3F4F6" stroke-width="0.7"/>',
        '</pattern>',
        '</defs>',
        f'<line x1="{panel_x - 24:.2f}" y1="68" x2="{panel_x - 24:.2f}" y2="{y_max + 10:.2f}" stroke="#D1D5DB" stroke-width="0.8"/>',
        text(panel_x, 34, "Gene-island neighborhoods", 15, "#111827", "700"),
        text(panel_x, 51, "Centered on the focal AioA/IdrA-related gene", 9, "#4B5563"),
    ]

    legend = [
        ("A", GENE_COLORS["AioA/IdrA-related A"]),
        ("B", GENE_COLORS["IdrB/AioB-related small subunit"]),
        ("P-like", GENE_COLORS["P-like"]),
        ("other CDS", GENE_COLORS["other annotated CDS"]),
        ("unknown", GENE_COLORS["function_unknown"]),
    ]
    lx = panel_x
    for label, color in legend:
        parts.append(
            f'<rect x="{lx:.2f}" y="61" width="14" height="8" fill="{color}" '
            'stroke="#374151" stroke-width="0.35"/>'
        )
        parts.append(text(lx + 19, 69, label, 8, "#374151"))
        lx += 96 if label not in {"other CDS", "unknown"} else 124

    scale_width = 5000.0 / (2 * extent) * panel_width
    scale_x = panel_x + panel_width - scale_width
    parts.extend(
        [
            f'<line x1="{scale_x:.2f}" y1="86" x2="{panel_x + panel_width:.2f}" y2="86" stroke="#111827" stroke-width="1.2"/>',
            f'<line x1="{scale_x:.2f}" y1="82" x2="{scale_x:.2f}" y2="90" stroke="#111827" stroke-width="1"/>',
            f'<line x1="{panel_x + panel_width:.2f}" y1="82" x2="{panel_x + panel_width:.2f}" y2="90" stroke="#111827" stroke-width="1"/>',
            text(scale_x + scale_width / 2, 80, "5 kb", 8, "#111827", "400", "middle"),
            f'<line x1="{center_x:.2f}" y1="{y_min - 9:.2f}" x2="{center_x:.2f}" y2="{y_max + 9:.2f}" stroke="#94A3B8" stroke-width="0.65" stroke-dasharray="3 3"/>',
        ]
    )

    for label in displayed:
        y = positions[label][1] - 2.2
        parts.append(
            f'<line x1="{panel_x:.2f}" y1="{y:.2f}" x2="{panel_x + panel_width:.2f}" y2="{y:.2f}" '
            'stroke="#E5E7EB" stroke-width="0.35"/>'
        )
        for row in rows_by_label[label]:
            try:
                x1 = sx(row["relative_bp_start"])
                x2 = sx(row["relative_bp_end"])
            except (KeyError, TypeError, ValueError):
                continue
            plot_label = row.get("plot_label", "other annotated CDS")
            fill = GENE_COLORS.get(plot_label, GENE_COLORS["other annotated CDS"])
            stroke = "#334155" if plot_label != "other annotated CDS" else "#94A3B8"
            parts.append(
                f'<polygon data-candidate-id="{esc(label)}" points="{arrow_points(x1, x2, y, 9, row.get("strand", "+"))}" '
                f'fill="{fill}" stroke="{stroke}" stroke-width="0.45"/>'
            )
            short = SHORT_LABELS.get(plot_label)
            if short and abs(x2 - x1) >= 11:
                parts.append(text((x1 + x2) / 2, y + 2.1, short, 5.8, "#FFFFFF", "700", "middle"))

    parts.append("</g>")
    return "\n".join(parts), displayed, extent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-svg", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-svg", type=Path, required=True)
    args = parser.parse_args()

    svg = args.input_svg.read_text()
    manifest = read_tsv(args.manifest)
    rows_by_label = defaultdict(list)
    for row in manifest:
        rows_by_label[row["tree_label"]].append(row)
    augment_r2111(rows_by_label)

    positions, source_labels = parse_tip_positions(svg, rows_by_label)
    missing = sorted(set(rows_by_label) - set(positions))
    if missing:
        raise SystemExit("Tree labels without SVG positions: " + ", ".join(missing))

    panel_x = 1190.0
    panel_width = 1010.0
    overlay, displayed, extent = build_overlay(rows_by_label, positions, panel_x, panel_width)
    svg = re.sub(
        r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"',
        lambda match: f'viewBox="0 0 2250 {match.group(2)}"',
        svg,
        count=1,
    )
    svg = svg.replace("</svg>", overlay + "\n</svg>")
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(svg)

    print(f"matched_gene_island_rows={len(displayed)}")
    print(f"relative_coordinate_extent=-{int(extent)}..{int(extent)} bp")
    print("manual_label_aliases=" + ", ".join(
        f"{source_labels[label]} -> {label}"
        for label in displayed
        if source_labels[label] != label
    ))
    print(f"output={args.output_svg}")


if __name__ == "__main__":
    main()
