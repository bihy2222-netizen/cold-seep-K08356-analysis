#!/usr/bin/env python3
import csv
import math
import re
from pathlib import Path


BASE = Path(__file__).resolve().parent
TREE = BASE / "AioA_IdrA_Unknown_group_derep49_plus_refs136.trimmed.fasta.treefile"
HITDATA = BASE / "hitdata.txt"
META = BASE / "beautified_small_tree" / "tip_metadata.tsv"
OLD_MANIFEST = BASE.parent / "MAG4-smallsmall" / "Fig5_style_tree_neighborhood_AoR_rooted" / "arrow_plot_gene_coordinate_manifest_48.tsv"
OUT_DIR = BASE / "small_tree_arrow_plot"

HABITAT_COLORS = {
    "IS": "#0072B2",
    "AS": "#E69F00",
    "ES": "#009E73",
    "NS": "#D55E00",
    "Reference": "#9CA3AF",
    "Unknown_ref": "#7B2CBF",
}

GENE_COLORS = {
    "AioA/IdrA-related A": "#2166AC",
    "IdrB/AioB-related small subunit": "#F58518",
    "canonical AioB": "#B279A2",
    "P-like": "#2CA25F",
    "function_unknown": "#6B7280",
    "other annotated CDS": "#D7DCE2",
    "CD-search only A": "#1F78B4",
}

KEY = {
    "NS|R2111_N500_0-10_bin13-k141_4404301_3": ("New NS K08356", "#D62728"),
    "IS|SY368YW-8-12_bin16-k141_314641_4": ("Sister IS K08356", "#FF9F1C"),
}


class Node:
    def __init__(self, name="", length=0.0, children=None):
        self.name = name
        self.length = length
        self.children = children or []
        self.x = 0.0
        self.y = 0.0

    @property
    def is_leaf(self):
        return not self.children


def esc(text):
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def parse_newick(text):
    text = text.strip().rstrip(";")
    i = 0

    def parse_name():
        nonlocal i
        start = i
        while i < len(text) and text[i] not in ":,()":
            i += 1
        return text[start:i].strip()

    def parse_length():
        nonlocal i
        if i >= len(text) or text[i] != ":":
            return 0.0
        i += 1
        start = i
        while i < len(text) and text[i] not in ",()":
            i += 1
        raw = text[start:i]
        try:
            return float(raw)
        except ValueError:
            return 0.0

    def parse_subtree():
        nonlocal i
        if text[i] == "(":
            i += 1
            children = []
            while True:
                children.append(parse_subtree())
                if i >= len(text):
                    break
                if text[i] == ",":
                    i += 1
                    continue
                if text[i] == ")":
                    i += 1
                    break
            name = parse_name()
            length = parse_length()
            return Node(name=name, length=length, children=children)
        name = parse_name()
        length = parse_length()
        return Node(name=name, length=length)

    return parse_subtree()


def walk(node):
    yield node
    for child in node.children:
        yield from walk(child)


def leaves(node):
    return [n for n in walk(node) if n.is_leaf]


def assign_tree_coords(root, y_by_label):
    def rec(node, x):
        node.x = x
        if node.is_leaf:
            node.y = y_by_label[node.name]
        else:
            for child in node.children:
                rec(child, x + child.length)
            node.y = sum(child.y for child in node.children) / len(node.children)

    rec(root, 0.0)


def strip_prefix(label):
    if label.startswith("ISMEJFig5|"):
        return label
    return label.split("|", 1)[1] if "|" in label else label


def habitat_for(label, meta):
    if label in meta and meta[label].get("is_mag") == "TRUE":
        return label.split("|", 1)[0]
    if label.startswith("ISMEJFig5|Unknown_clade|"):
        return "Unknown_ref"
    return "Reference"


def read_meta():
    if not META.exists():
        return {}
    with META.open(newline="") as handle:
        return {r["label"]: r for r in csv.DictReader(handle, delimiter="\t")}


def read_manifest():
    rows_by_id = {}
    if not OLD_MANIFEST.exists():
        return rows_by_id
    with OLD_MANIFEST.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows_by_id.setdefault(row["candidate_id"], []).append(row)
    return rows_by_id


def read_hitdata_focus():
    focus = {}
    with HITDATA.open() as handle:
        for line in handle:
            if not line.startswith("Q#"):
                continue
            parts = line.rstrip("\n").split("\t")
            q = parts[0]
            m = re.search(r">\s*(.+?)\s+#\s*(\d+)\s+#\s*(\d+)\s+#\s*(-?1)", q)
            if not m:
                continue
            label, start, end, strand = m.groups()
            acc = parts[7] if len(parts) > 7 else ""
            short = parts[8] if len(parts) > 8 else ""
            if label not in focus or short == "arsenite_ox_L superfamily":
                focus[label] = {
                    "label": label,
                    "candidate_id": strip_prefix(label),
                    "start": int(start),
                    "end": int(end),
                    "strand": "+" if strand == "1" else "-",
                    "annotation": short or acc or "CD-search hit",
                }
    return focus


def gene_label(row):
    family = row.get("joint_reference_family", "none")
    plot = row.get("plot_label", "")
    if family == "AioA/IdrA-related A":
        return "AioA/IdrA-related A"
    if family == "IdrB_related":
        return "IdrB/AioB-related small subunit"
    if family == "canonical_AioB":
        return "canonical AioB"
    if family == "P_like":
        return "P-like"
    if family == "function_unknown":
        return "function_unknown"
    if plot == "AioA/IdrA-related A":
        return "AioA/IdrA-related A"
    return "other annotated CDS"


def arrow_path(x1, x2, y, h, strand):
    x1, x2 = float(x1), float(x2)
    if x2 < x1:
        x1, x2 = x2, x1
    width = max(x2 - x1, 5.0)
    head = min(10.0, width * 0.42)
    y1, y2, ym = y - h / 2, y + h / 2, y
    if strand == "-":
        pts = [(x1, ym), (x1 + head, y1), (x2, y1), (x2, y2), (x1 + head, y2)]
    else:
        pts = [(x1, y1), (x2 - head, y1), (x2, ym), (x2 - head, y2), (x1, y2)]
    return " ".join(f"{a:.1f},{b:.1f}" for a, b in pts)


def write_svg(root, meta, manifest, focus, mag_only=False, filename="small_tree_with_gene_arrows.svg"):
    all_leaves = leaves(root)
    if mag_only:
        all_leaves = [leaf for leaf in all_leaves if leaf.name in focus]
    row_h = 15
    top = 110
    leaf_y = {leaf.name: top + i * row_h for i, leaf in enumerate(all_leaves)}
    if mag_only:
        for leaf in all_leaves:
            leaf.y = leaf_y[leaf.name]
    else:
        assign_tree_coords(root, leaf_y)

    max_x = max(n.x for n in walk(root)) or 1.0
    left = 32
    label_x = 74
    arrow_x = 560
    arrow_w = 760
    evidence_x = arrow_x + arrow_w + 26
    width = evidence_x + 290
    height = top + len(all_leaves) * row_h + 92

    mag_labels = [leaf.name for leaf in all_leaves if leaf.name in focus]
    rel_min, rel_max = -12000, 14500
    for label in mag_labels:
        rows = manifest.get(strip_prefix(label), [])
        if rows:
            vals = []
            for r in rows:
                try:
                    vals.extend([int(r["relative_bp_start"]), int(r["relative_bp_end"])])
                except Exception:
                    pass
            if vals:
                rel_min = min(rel_min, min(vals))
                rel_max = max(rel_max, max(vals))

    def ax(rel):
        return arrow_x + (rel - rel_min) / (rel_max - rel_min) * arrow_w

    chunks = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        '<text x="32" y="30" font-family="Arial" font-size="20" font-weight="700" fill="#111827">Small AioA/IdrA MAG gene-arrow annotations</text>',
        '<text x="32" y="52" font-family="Arial" font-size="12" fill="#4B5563">Rows follow the iTOL/treefile top-to-bottom tip order approximately; tree topology is intentionally hidden. Color strips mark habitats.</text>',
        '<text x="32" y="82" font-family="Arial" font-size="11" font-weight="700" fill="#111827">Sequence</text>',
        f'<text x="{arrow_x}" y="82" font-family="Arial" font-size="11" font-weight="700" fill="#111827">Gene neighborhood / CD-search arrows</text>',
        f'<text x="{evidence_x}" y="82" font-family="Arial" font-size="11" font-weight="700" fill="#111827">Evidence</text>',
    ]

    lx = 32
    for label, color in [
        ("IS", HABITAT_COLORS["IS"]),
        ("AS", HABITAT_COLORS["AS"]),
        ("ES", HABITAT_COLORS["ES"]),
        ("NS", HABITAT_COLORS["NS"]),
        ("Unknown ref", HABITAT_COLORS["Unknown_ref"]),
        ("Reference", HABITAT_COLORS["Reference"]),
    ]:
        chunks.append(f'<rect x="{lx}" y="92" width="12" height="8" fill="{color}"/>')
        chunks.append(f'<text x="{lx+16}" y="100" font-family="Arial" font-size="10" fill="#374151">{esc(label)}</text>')
        lx += 88 if label != "Unknown ref" else 118
    for label, color in [
        ("AioA/IdrA A", GENE_COLORS["AioA/IdrA-related A"]),
        ("IdrB/AioB", GENE_COLORS["IdrB/AioB-related small subunit"]),
        ("P-like", GENE_COLORS["P-like"]),
        ("canonical AioB", GENE_COLORS["canonical AioB"]),
        ("other CDS", GENE_COLORS["other annotated CDS"]),
    ]:
        chunks.append(f'<rect x="{lx}" y="92" width="12" height="8" fill="{color}" stroke="#374151" stroke-width="0.4"/>')
        chunks.append(f'<text x="{lx+16}" y="100" font-family="Arial" font-size="10" fill="#374151">{esc(label)}</text>')
        lx += 112

    chunks.append(f'<line x1="{ax(0):.1f}" y1="{top-10}" x2="{ax(0):.1f}" y2="{height-72}" stroke="#94A3B8" stroke-width="0.8" stroke-dasharray="3 4"/>')
    chunks.append(f'<text x="{ax(0):.1f}" y="{height-52}" text-anchor="middle" font-family="Arial" font-size="10" fill="#64748B">focal AioA/IdrA</text>')

    combined_rows = []
    for leaf in all_leaves:
        label = leaf.name
        y = leaf.y
        habitat = habitat_for(label, meta)
        color = KEY.get(label, (None, None))[1] or HABITAT_COLORS.get(habitat, "#9CA3AF")
        size = 10 if label in focus or label.startswith("ISMEJFig5|Unknown_clade|") else 7
        weight = "700" if label in KEY or label.startswith("ISMEJFig5|Unknown_clade|") else "400"
        display = meta.get(label, {}).get("display_label", label)
        if len(display) > 74:
            display = display[:71] + "..."
        chunks.append(f'<rect x="{left}" y="{y-5:.2f}" width="28" height="10" rx="1.5" fill="{color}" opacity="0.95"/>')
        chunks.append(f'<text x="{label_x}" y="{y+3:.2f}" font-family="Arial" font-size="{size}" font-weight="{weight}" fill="{color}">{esc(display)}</text>')

        if label not in focus:
            continue
        cid = strip_prefix(label)
        rows = manifest.get(cid, [])
        chunks.append(f'<rect x="{arrow_x-8}" y="{y-6.5:.1f}" width="5" height="13" rx="1" fill="{HABITAT_COLORS.get(habitat, "#9CA3AF")}"/>')
        if rows:
            for r in rows:
                try:
                    x1 = ax(int(r["relative_bp_start"]))
                    x2 = ax(int(r["relative_bp_end"]))
                except Exception:
                    continue
                glabel = gene_label(r)
                fill = GENE_COLORS[glabel]
                stroke = "#111827" if glabel != "other annotated CDS" else "#94A3B8"
                chunks.append(f'<polygon points="{arrow_path(x1, x2, y, 8.8, r.get("strand", "+"))}" fill="{fill}" stroke="{stroke}" stroke-width="0.45"/>')
                combined_rows.append({
                    "tree_label": label,
                    "candidate_id": cid,
                    "source": r.get("evidence_source", "previous manifest"),
                    "gene_id": r.get("gene_id", ""),
                    "relative_gene_position": r.get("relative_gene_position", ""),
                    "relative_bp_start": r.get("relative_bp_start", ""),
                    "relative_bp_end": r.get("relative_bp_end", ""),
                    "strand": r.get("strand", ""),
                    "plot_label": glabel,
                    "neighborhood_structure_class": r.get("neighborhood_structure_class", ""),
                    "strict_synteny_QC_pass": r.get("strict_synteny_QC_pass", ""),
                    "annotation": r.get("DRAM_annotation", ""),
                })
            focal = next((r for r in rows if r.get("relative_gene_position") == "0"), rows[0])
            evidence = focal.get("neighborhood_structure_class", "")
            qc = focal.get("strict_synteny_QC_pass", "")
            chunks.append(f'<text x="{evidence_x}" y="{y+3:.2f}" font-family="Arial" font-size="9" fill="#374151">{esc(evidence)}; QC={esc(qc)}</text>')
        else:
            f = focus[label]
            rel_half = max((f["end"] - f["start"]) / 2, 200)
            chunks.append(f'<polygon points="{arrow_path(ax(-rel_half), ax(rel_half), y, 9.5, f["strand"])}" fill="{GENE_COLORS["CD-search only A"]}" stroke="#111827" stroke-width="0.55"/>')
            chunks.append(f'<text x="{evidence_x}" y="{y+3:.2f}" font-family="Arial" font-size="9" fill="#374151">CD-search only; no neighborhood manifest</text>')
            combined_rows.append({
                "tree_label": label,
                "candidate_id": cid,
                "source": "hitdata CD-search only",
                "gene_id": cid,
                "relative_gene_position": "0",
                "relative_bp_start": int(-rel_half),
                "relative_bp_end": int(rel_half),
                "strand": f["strand"],
                "plot_label": "CD-search only A",
                "neighborhood_structure_class": "not available",
                "strict_synteny_QC_pass": "not available",
                "annotation": f["annotation"],
            })
        if label in KEY:
            key_label, key_color = KEY[label]
            chunks.append(f'<text x="{arrow_x-148}" y="{y+3:.2f}" font-family="Arial" font-size="9" font-weight="700" fill="{key_color}">{esc(key_label)}</text>')

    chunks.append(f'<text x="32" y="{height-24}" font-family="Arial" font-size="10" fill="#64748B">Note: 48 MAG neighborhoods are drawn from the previous coordinate manifest; R2111_N500_0-10_bin13 is drawn from CD-search focus coordinates only because no neighborhood row was present in that manifest.</text>')
    chunks.append("</svg>")

    OUT_DIR.mkdir(exist_ok=True)
    svg_path = OUT_DIR / filename
    svg_path.write_text("\n".join(chunks))

    tsv_path = OUT_DIR / "small_tree_gene_arrow_manifest.tsv"
    with tsv_path.open("w", newline="") as handle:
        fields = [
            "tree_label",
            "candidate_id",
            "source",
            "gene_id",
            "relative_gene_position",
            "relative_bp_start",
            "relative_bp_end",
            "strand",
            "plot_label",
            "neighborhood_structure_class",
            "strict_synteny_QC_pass",
            "annotation",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(combined_rows)

    readme = OUT_DIR / "README_small_tree_arrow_plot.txt"
    readme.write_text(
        "Small tree arrow plot outputs\n"
        "=============================\n"
        f"Treefile: {TREE}\n"
        f"CD-search hitdata: {HITDATA}\n"
        f"Previous coordinate manifest reused: {OLD_MANIFEST}\n\n"
        "Main figure:\n"
        f"- {svg_path}\n\n"
        "Data table:\n"
        f"- {tsv_path}\n\n"
        "Important note: the previous neighborhood coordinate manifest covered 48 MAG candidates. "
        "The new NS candidate R2111_N500_0-10_bin13-k141_4404301_3 was present in hitdata but absent "
        "from that manifest, so only its focal CD-search AioA/IdrA arrow is drawn.\n"
    )
    return svg_path, tsv_path, readme


def main():
    root = parse_newick(TREE.read_text())
    meta = read_meta()
    manifest = read_manifest()
    focus = read_hitdata_focus()
    svg, tsv, readme = write_svg(root, meta, manifest, focus, mag_only=False, filename="small_tree_all_tips_gene_arrows_iTOL_order.svg")
    mag_svg, _, _ = write_svg(root, meta, manifest, focus, mag_only=True, filename="small_tree_MAG_only_gene_arrows_iTOL_order.svg")
    print(mag_svg)
    print(svg)
    print(tsv)
    print(readme)


if __name__ == "__main__":
    main()
