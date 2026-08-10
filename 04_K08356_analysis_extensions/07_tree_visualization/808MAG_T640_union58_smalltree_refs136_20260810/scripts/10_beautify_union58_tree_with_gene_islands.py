#!/usr/bin/env python3
"""Beautify the supplied latest tree SVG and add aligned 58-candidate gene islands."""

import csv
import html
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
INPUT_SVG = ROOT / "00_input/latest_union58_itol_export.svg"
EVIDENCE = ROOT / "04_tables/candidate58_four_clade_integrated_evidence.tsv"
MANIFEST = ROOT / "04_neighborhood/union58_gene_island_manifest.tsv"
OUTDIR = ROOT / "figures"
OUTPUT_SVG = OUTDIR / "MAG_T640_union58_four_clades_gene_islands.svg"

CLASS_ORDER = ["Canonical AioA clade", "Unknown DMSOR clade", "IdrA-associated clade", "IdrA clade"]
CLASS_COLORS = {
    "Canonical AioA clade": "#009E73",
    "Unknown DMSOR clade": "#D55E00",
    "IdrA-associated clade": "#0072B2",
    "IdrA clade": "#7B3294",
}
GENE_COLORS = {
    "Molybdopterin oxidoreductase": "#C92C2C",
    "Rieske [2Fe-2S] small subunit": "#D89B18",
    "IdrP-like accessory protein": "#3A9D4E",
    "Cytochrome c / peroxidase": "#3E4F8A",
    "Redox / oxidoreductase": "#A65A2E",
    "Transporter": "#8E4E9B",
    "Arsenic resistance protein": "#B87838",
    "Sulfur oxidation protein": "#2A8F8A",
    "Other annotated protein": "#D7DCE2",
    "Non-conserved protein": "#E7E7E4",
}


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def esc(value):
    return html.escape(str(value), quote=True)


def aliases(label):
    candidate = label.split("|", 1)[1] if "|" in label else label
    return {
        label, candidate,
        label.replace("_", " "), candidate.replace("_", " "),
        label.replace("_bin", " bin"), candidate.replace("_bin", " bin"),
    }


def arrow_points(x1, x2, y, height, strand):
    if x2 < x1:
        x1, x2 = x2, x1
    width = max(x2 - x1, 5)
    head = min(13, width * 0.38)
    top, bottom = y - height / 2, y + height / 2
    if strand == "-":
        pts = [(x1, y), (x1 + head, top), (x2, top), (x2, bottom), (x1 + head, bottom)]
    else:
        pts = [(x1, top), (x2 - head, top), (x2, y), (x2 - head, bottom), (x1, bottom)]
    return " ".join(f"{x:.1f},{yy:.1f}" for x, yy in pts)


evidence = read_tsv(EVIDENCE)
manifest = read_tsv(MANIFEST)
by_tree_label = {row["prefixed_target_ID"]: row for row in evidence}
alias_map = {}
for label in by_tree_label:
    for alias in aliases(label):
        alias_map[alias] = label

svg = INPUT_SVG.read_text()
root_end = svg.find(">")
svg = svg[:root_end + 1] + '<rect id="white-background" width="1600" height="1180" fill="#ffffff"/>\n' + svg[root_end + 1:]
svg = re.sub(r'width="1359"', 'width="1600"', svg, count=1)
svg = re.sub(r'height="724"', 'height="1180"', svg, count=1)
svg = re.sub(r'viewBox="0,0,1359,724"', 'viewBox="0,0,1600,1180"', svg, count=1)

text_pattern = re.compile(r'<text\b(?P<attrs>[^>]*\bx="(?P<x>[-0-9.]+)"[^>]*\by="(?P<y>[-0-9.]+)"[^>]*)>(?P<label>.*?)</text>', re.S)
positions = {}

def color_tip(match):
    visible = re.sub(r"<.*?>", "", match.group("label")).strip()
    canonical = alias_map.get(visible)
    if not canonical:
        return match.group(0)
    row = by_tree_label[canonical]
    positions[canonical] = (float(match.group("x")), float(match.group("y")))
    attrs = re.sub(r'\sfill="[^"]*"', "", match.group("attrs"))
    attrs = re.sub(r'\sfont-weight="[^"]*"', "", attrs)
    attrs += f' fill="{CLASS_COLORS[row["Final_class"]]}" font-weight="bold" data-final-class="{esc(row["Final_class"])}"'
    return f'<text{attrs}>{match.group("label")}</text>'

svg = text_pattern.sub(color_tip, svg)
missing = sorted(set(by_tree_label) - set(positions))
if missing:
    raise SystemExit("Candidate labels not found in supplied SVG: " + ", ".join(missing))

rows_by_label = defaultdict(list)
for row in manifest:
    rows_by_label[row["tree_label"]].append(row)

extent = max(abs(float(row[key])) for row in manifest for key in ("relative_bp_start", "relative_bp_end"))
extent = max(10000, math.ceil(extent / 1000) * 1000)
panel_x, panel_width = 1775.0, 2100.0
center_x = panel_x + panel_width / 2

def sx(value):
    return center_x + float(value) / (2 * extent) * panel_width

y_values = [pos[1] for pos in positions.values()]
parts = [
    '<g id="union58-gene-island-panel" font-family="Times New Roman, Times, serif">',
    f'<line x1="{center_x:.1f}" y1="{min(y_values)-12:.1f}" x2="{center_x:.1f}" y2="{max(y_values)+12:.1f}" stroke="#94A3B8" stroke-width="1.2" stroke-dasharray="6,5"/>',
]

for label, (x, y) in sorted(positions.items(), key=lambda item: item[1][1]):
    cls = by_tree_label[label]["Final_class"]
    parts.append(f'<rect x="1625" y="{y-10:.1f}" width="20" height="18" rx="2" fill="{CLASS_COLORS[cls]}" stroke="#ffffff" stroke-width="1.2"/>')
    parts.append(f'<line x1="1660" y1="{y-1:.1f}" x2="{panel_x+panel_width:.1f}" y2="{y-1:.1f}" stroke="#E5E7EB" stroke-width="0.55"/>')
    for gene in rows_by_label[label]:
        x1, x2 = sx(gene["relative_bp_start"]), sx(gene["relative_bp_end"])
        plot_label = gene["plot_label"]
        fill = GENE_COLORS[plot_label]
        stroke = "#4B5563" if plot_label != "Non-conserved protein" else "#777777"
        parts.append(f'<polygon data-candidate="{esc(label)}" data-gene="{esc(gene["gene_id"])}" points="{arrow_points(x1,x2,y-1,15,gene["strand"])}" fill="{fill}" stroke="{stroke}" stroke-width="0.7"/>')

scale_width = 5000 / (2 * extent) * panel_width
parts.append('</g>')
gene_overlay = "\n".join(parts)
needle = '</g><g id="wrapperHolder"'
if needle not in svg:
    raise SystemExit("Could not locate end of mainHolder in SVG")
svg = svg.replace(needle, gene_overlay + '\n</g><g id="wrapperHolder"', 1)

counts = Counter(row["Final_class"] for row in evidence)
legend = ['<g id="four-clade-and-gene-legends" font-family="Times New Roman, Times, serif">',
          '<text x="360" y="18" font-size="14" font-weight="bold" fill="#111827">Candidate classification</text>']
x = 360
for cls in CLASS_ORDER:
    legend.append(f'<rect x="{x}" y="27" width="12" height="10" fill="{CLASS_COLORS[cls]}"/>')
    short_cls = {"Canonical AioA clade":"Canonical AioA", "Unknown DMSOR clade":"Unknown DMSOR", "IdrA-associated clade":"IdrA-associated", "IdrA clade":"IdrA"}[cls]
    legend.append(f'<text x="{x+17}" y="36" font-size="10" fill="#1F2937">{esc(short_cls)} (n={counts[cls]})</text>')
    x += 155
legend.append('<text x="1010" y="13" font-size="14" font-weight="bold" fill="#111827">Protein subfamily</text>')
legend_order = list(GENE_COLORS)
for i, label in enumerate(legend_order):
    col, row = i // 5, i % 5
    x, y = 1010 + col * 285, 22 + row * 14
    legend.append(f'<polygon points="{x},{y} {x+9},{y+5} {x},{y+10}" fill="{GENE_COLORS[label]}" stroke="#4B5563" stroke-width="0.45"/>')
    legend.append(f'<text x="{x+14}" y="{y+8}" font-size="9.5" fill="#1F2937">{esc(label)}</text>')
legend.extend([
    '<line x1="1510" y1="70" x2="1570" y2="70" stroke="#111827" stroke-width="1.2"/>',
    '<line x1="1510" y1="66" x2="1510" y2="74" stroke="#111827" stroke-width="1"/>',
    '<line x1="1570" y1="66" x2="1570" y2="74" stroke="#111827" stroke-width="1"/>',
    '<text x="1540" y="63" text-anchor="middle" font-size="9" fill="#111827">5 kb</text>',
])
legend.append('</g>')
svg = svg.replace('</svg>', "\n".join(legend) + '\n</svg>')

OUTDIR.mkdir(parents=True, exist_ok=True)
OUTPUT_SVG.write_text(svg)
print(f"candidate_positions={len(positions)}")
print(f"gene_island_rows={len(manifest)}")
print(f"coordinate_extent=+/-{int(extent)} bp")
print(f"output={OUTPUT_SVG}")
