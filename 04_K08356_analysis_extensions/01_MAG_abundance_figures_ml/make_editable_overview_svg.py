#!/usr/bin/env python3
import csv
import html
import os
import re
from collections import OrderedDict

BASE = "/Users/catherine/Downloads/MAG-bin-tpm/K08356_MAG_gene_island_arrows_20260708"
COORDS = os.path.join(BASE, "K08356_four_branch_gene_neighborhood_coordinates.csv")
FIG_SUBDIR = os.environ.get("K08356_FIG_SUBDIR", "four_branch_gene_island_figures")
SHOW_DOMAIN_MARKERS = os.environ.get("K08356_SHOW_DOMAIN_MARKERS", "true").lower() in {"1", "true", "yes", "y"}
SHOW_GENE_TEXT_LABELS = os.environ.get("K08356_SHOW_GENE_TEXT_LABELS", "true").lower() in {"1", "true", "yes", "y"}
SELECTED_IDS_FILE = os.environ.get("K08356_SELECTED_IDS_FILE", "")
OUT = os.path.join(
    BASE,
    FIG_SUBDIR,
    "K08356_four_branch_MAG_gene_island_arrows_all_editable.svg",
)
os.makedirs(os.path.dirname(OUT), exist_ok=True)

WIDTH, HEIGHT = 1180, 1550
LEFT, RIGHT = 260, 920
TOP, BOTTOM = 92, 1370
STRIP_X0, STRIP_X1 = 930, 956
LEGEND_X = 975
FONT = "Times New Roman, Times, serif"

BRANCH_ORDER = [
    "IdrA-associated",
    "canonical aioA-associated",
    "aioA-like-associated",
    "unknown AioA-like / uncertain DMSOR",
]

COLORS = {
    "Target IdrA-like K08356": "#5B8DB8",
    "Target canonical aioA": "#2CA25F",
    "Target aioA-like K08356": "#009E73",
    "Target uncertain AioA-like": "#CC79A7",
    "Rieske/AioB": "#B9A7D8",
    "Cytochrome/peroxidase": "#9E9E9E",
    "Arsenic resistance": "#A8A8A8",
    "Transporter": "#BDBDBD",
    "Nitrogen metabolism": "#969696",
    "Sox/sulfur oxidation": "#8A8A8A",
    "Redox/oxidoreductase": "#C7C7C7",
    "Annotated other": "#D6D6D6",
    "Hypothetical/unknown": "#EFEFEF",
}

DOMAIN_MARKERS = {
    "target": ("Molybdopterin oxidoreductase / K08356", "#E41A1C"),
    "rieske": ("Rieske [2Fe-2S] / AioB", "#FFB000"),
    "cyt": ("Cytochrome c peroxidase", "#009E3D"),
}


def esc(x):
    return html.escape(str(x), quote=True)


def arrow_points(x1, x2, y, direction, h=1.9):
    length = max(x2 - x1, 1)
    head = min(length * 0.18, length * 0.45, 10)
    if int(float(direction)) >= 0:
        pts = [
            (x1, y - h),
            (x2 - head, y - h),
            (x2 - head, y - h * 1.55),
            (x2, y),
            (x2 - head, y + h * 1.55),
            (x2 - head, y + h),
            (x1, y + h),
        ]
    else:
        pts = [
            (x2, y - h),
            (x1 + head, y - h),
            (x1 + head, y - h * 1.55),
            (x1, y),
            (x1 + head, y + h * 1.55),
            (x1 + head, y + h),
            (x2, y + h),
        ]
    return " ".join(f"{x:.2f},{yy:.2f}" for x, yy in pts)


rows = []
selected_ids = None
if SELECTED_IDS_FILE:
    with open(SELECTED_IDS_FILE) as handle:
        selected_ids = {line.strip() for line in handle if line.strip()}

with open(COORDS, newline="") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        if selected_ids is not None and row["Original_MAG_ID"] not in selected_ids:
            continue
        row["Start"] = float(row["Start"])
        row["End"] = float(row["End"])
        row["Direction"] = int(float(row["Direction"]))
        row["Molecule_short"] = f"{row['Original_MAG_ID']} | {row['Corrected_habitat']}"
        rows.append(row)

for mol in {r["Molecule_short"] for r in rows}:
    mol_rows = [r for r in rows if r["Molecule_short"] == mol]
    region_start = min(r["Start"] for r in mol_rows)
    for r in mol_rows:
        r["Start_plot"] = r["Start"] - region_start + 1
        r["End_plot"] = r["End"] - region_start + 1
        r["Mid_plot"] = (r["Start_plot"] + r["End_plot"]) / 2

x_max = max(r["End_plot"] for r in rows)

by_branch = OrderedDict((b, []) for b in BRANCH_ORDER)
for b in BRANCH_ORDER:
    seen = OrderedDict()
    for r in rows:
        if r["Feature"] == b and r["Molecule_short"] not in seen:
            seen[r["Molecule_short"]] = None
    by_branch[b] = list(seen.keys())

total_rows = sum(len(v) for v in by_branch.values())
gap = 20
row_step = (BOTTOM - TOP - gap * (len(BRANCH_ORDER) - 1)) / max(total_rows - 1, 1)

y_by_mol = {}
branch_bounds = {}
y = TOP
for b in BRANCH_ORDER:
    start_y = y - row_step * 0.55
    for mol in by_branch[b]:
        y_by_mol[mol] = y
        y += row_step
    end_y = y - row_step * 0.45
    branch_bounds[b] = (start_y, end_y)
    y += gap

def sx(x):
    return LEFT + (x / x_max) * (RIGHT - LEFT)

parts = [
    '<?xml version="1.0" encoding="UTF-8"?>',
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}pt" height="{HEIGHT}pt" viewBox="0 0 {WIDTH} {HEIGHT}">',
    '<rect width="100%" height="100%" fill="white"/>',
    f'<text x="{WIDTH/2}" y="26" text-anchor="middle" font-family="{FONT}" font-size="18" font-weight="bold">K08356 MAG gene neighborhoods across four protein-branch groups</text>',
    f'<text x="{WIDTH/2}" y="45" text-anchor="middle" font-family="{FONT}" font-size="10">Target-centered neighborhoods from downloaded sample-level CDS files; up to 10 upstream and 10 downstream ORFs are shown</text>',
]

if SHOW_DOMAIN_MARKERS:
    marker_y = 76
    parts.append(f'<text x="{LEGEND_X}" y="{marker_y}" font-family="{FONT}" font-size="9" font-weight="bold">Functional domain markers</text>')
    for idx, (_, (lab, col)) in enumerate(DOMAIN_MARKERS.items()):
        yy = marker_y + 18 + idx * 18
        parts.append(f'<circle cx="{LEGEND_X+6}" cy="{yy-3}" r="3.2" fill="{col}" stroke="black" stroke-width="0.5"/>')
        parts.append(f'<text x="{LEGEND_X+18}" y="{yy}" font-family="{FONT}" font-size="7.6">{esc(lab)}</text>')

# grid and axis
for tick in range(0, int(x_max // 10000) * 10000 + 1, 10000):
    x = sx(tick)
    parts.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{TOP-10}" y2="{BOTTOM+5}" stroke="#e8e8e8" stroke-dasharray="5,5" stroke-width="1"/>')
    parts.append(f'<text x="{x:.2f}" y="{BOTTOM+22}" text-anchor="middle" font-family="{FONT}" font-size="8">{tick//1000}k</text>')
parts.append(f'<line x1="{LEFT}" x2="{RIGHT}" y1="{BOTTOM+8}" y2="{BOTTOM+8}" stroke="black" stroke-width="1"/>')
parts.append(f'<text x="{(LEFT+RIGHT)/2}" y="{BOTTOM+43}" text-anchor="middle" font-family="{FONT}" font-size="10" font-weight="bold">Relative genomic position in target-centered neighborhood (bp)</text>')

# strips
for b, (y0, y1) in branch_bounds.items():
    parts.append(f'<rect x="{STRIP_X0}" y="{y0:.2f}" width="{STRIP_X1-STRIP_X0}" height="{y1-y0:.2f}" fill="#f2f2f2" stroke="#aaaaaa" stroke-width="0.6"/>')
    cx, cy = (STRIP_X0 + STRIP_X1) / 2, (y0 + y1) / 2
    parts.append(f'<text x="{cx:.2f}" y="{cy:.2f}" text-anchor="middle" font-family="{FONT}" font-size="8" font-weight="bold" transform="rotate(90 {cx:.2f} {cy:.2f})">{esc(b)}</text>')

# baselines
for mol, yv in y_by_mol.items():
    mol_rows = [r for r in rows if r["Molecule_short"] == mol]
    display = mol.replace(" | ", " | ")
    parts.append(f'<text x="{LEFT-8}" y="{yv+2.4:.2f}" text-anchor="end" font-family="{FONT}" font-size="6.3">{esc(display)}</text>')
    parts.append(f'<line x1="{sx(min(r["Start_plot"] for r in mol_rows)):.2f}" x2="{sx(max(r["End_plot"] for r in mol_rows)):.2f}" y1="{yv:.2f}" y2="{yv:.2f}" stroke="#bdbdbd" stroke-width="0.5"/>')

# arrows and labels
for i, r in enumerate(rows):
    yv = y_by_mol[r["Molecule_short"]]
    color = COLORS.get(r["Gene_class"], "#EFEFEF")
    x1, x2 = sx(r["Start_plot"]), sx(r["End_plot"])
    parts.append(
        f'<polygon id="orf_{i+1}" points="{arrow_points(x1, x2, yv, r["Direction"])}" '
        f'fill="{color}" stroke="#555555" stroke-width="0.35"/>'
    )
    marker_color = None
    if r.get("Is_target") == "yes":
        marker_color = DOMAIN_MARKERS["target"][1]
    elif r.get("Gene_class") == "Rieske/AioB":
        marker_color = DOMAIN_MARKERS["rieske"][1]
    elif r.get("Gene_class") == "Cytochrome/peroxidase":
        marker_color = DOMAIN_MARKERS["cyt"][1]
    if SHOW_DOMAIN_MARKERS and marker_color:
        parts.append(
            f'<circle cx="{sx(r["Mid_plot"]):.2f}" cy="{yv-4.8:.2f}" r="2.2" '
            f'fill="{marker_color}" stroke="black" stroke-width="0.45"/>'
        )
    label = r.get("Gene_label", "")
    if label and SHOW_GENE_TEXT_LABELS:
        label_y = yv - 9.5 if SHOW_DOMAIN_MARKERS else yv - 6.5
        parts.append(
            f'<text x="{sx(r["Mid_plot"]):.2f}" y="{label_y:.2f}" text-anchor="middle" '
            f'font-family="{FONT}" font-size="4.0">{esc(label)}</text>'
        )

# legend
legend_items = [(k, v) for k, v in COLORS.items()]
legend_y = 170
legend_x = LEGEND_X
parts.append(f'<text x="{legend_x}" y="{legend_y}" font-family="{FONT}" font-size="9" font-weight="bold">Gene/function class</text>')
for idx, (lab, col) in enumerate(legend_items):
    x = legend_x
    yy = legend_y + 16 + idx * 15
    parts.append(f'<rect x="{x}" y="{yy-10}" width="11" height="9" fill="{col}" stroke="#555555" stroke-width="0.35"/>')
    parts.append(f'<text x="{x+17}" y="{yy-2}" font-family="{FONT}" font-size="7.2">{esc(lab)}</text>')

caption_y = HEIGHT - 28
parts.append(f'<text x="{WIDTH/2}" y="{caption_y}" text-anchor="middle" font-family="{FONT}" font-size="8.5">Colored arrows indicate K08356 branch targets; grey arrows indicate neighboring ORFs. Circle markers denote key functional domains.</text>')

parts.append("</svg>")

with open(OUT, "w") as handle:
    handle.write("\n".join(parts))

print(OUT)
