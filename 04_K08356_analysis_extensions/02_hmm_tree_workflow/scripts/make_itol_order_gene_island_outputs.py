#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path
from xml.etree import ElementTree as ET


BASE = Path(__file__).resolve().parent
DEFAULT_INPUT_SVG = Path("/Users/catherine/Downloads/c0J-IDrhpB2odQ8jRmUo-Q.svg")
DEFAULT_TREEFILE = BASE / "AioA_IdrA_Unknown_group_derep49_plus_refs136.trimmed.fasta.treefile"
DEFAULT_OUT = BASE / "small_tree_itol_order_gene_island_outputs"
DEFAULT_META = BASE / "beautified_small_tree" / "tip_metadata.tsv"
DEFAULT_ARROW_MANIFEST = BASE / "small_tree_arrow_plot" / "small_tree_gene_arrow_manifest.tsv"
DEFAULT_HMM_SCORES = BASE.parent / "HMM_supplement_result_with_neighborhood_unpacked" / "result" / "03_tables" / "HMM_scores_48.tsv"
DEFAULT_T640 = BASE.parent / "HMM_supplement_result_with_neighborhood_unpacked" / "result" / "03_tables" / "all_MAG_T640_hits.tsv"
DEFAULT_STRICT_QC = BASE.parent / "joint_reference_search_48" / "strict_synteny_qc_48_reclassified.tsv"

HABITAT_COLORS = {
    "IS": "#0072B2",
    "AS": "#E69F00",
    "ES": "#009E73",
    "NS": "#D55E00",
    "Reference": "#8A8F98",
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


def esc(text):
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def read_tsv(path):
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


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


def parse_svg_tip_order(svg_text):
    rows = []
    for m in re.finditer(r"<text\b([^>]*)>(.*?)</text>", svg_text, re.S):
        attrs = m.group(1)
        value = re.sub(r"<.*?>", "", m.group(2)).strip()
        xm = re.search(r'\bx="([-0-9.]+)"', attrs)
        ym = re.search(r'\by="([-0-9.]+)"', attrs)
        if not (xm and ym and value):
            continue
        x = float(xm.group(1))
        y = float(ym.group(1))
        if x >= 900 and value not in {"out group"}:
            rows.append({"label": value, "x": x, "y": y})
    rows.sort(key=lambda r: r["y"])
    return rows


def parse_newick_leaf_order(treefile):
    text = treefile.read_text(errors="ignore").strip()
    return [m.group(1) for m in re.finditer(r"(?<=[(,])([^():,]+):", text)]


def augment_order_from_treefile(order, manifest_labels, treefile):
    if not treefile.exists():
        return order
    newick_labels = parse_newick_leaf_order(treefile)
    visible = {row["label"] for row in order}
    missing = [label for label in manifest_labels if label not in visible]
    if not missing:
        return order

    y_by_label = {row["label"]: float(row["y"]) for row in order}
    x_default = max([float(row["x"]) for row in order] + [1005.0])
    additions = []
    for label in missing:
        if label not in newick_labels:
            continue
        idx = newick_labels.index(label)
        prev_y = next((y_by_label[newick_labels[j]] for j in range(idx - 1, -1, -1) if newick_labels[j] in y_by_label), None)
        next_y = next((y_by_label[newick_labels[j]] for j in range(idx + 1, len(newick_labels)) if newick_labels[j] in y_by_label), None)
        if prev_y is not None and next_y is not None:
            y = prev_y + (next_y - prev_y) / 2
        elif prev_y is not None:
            y = prev_y + 30
        elif next_y is not None:
            y = next_y - 30
        else:
            y = max([float(row["y"]) for row in order] + [0.0]) + 30
        while any(abs(float(row["y"]) - y) < 0.01 for row in order + additions):
            y += 0.1
        additions.append({"label": label, "x": x_default, "y": y, "order_note": "inserted_from_newick_neighbor_order"})
    combined = [dict(row, order_note=row.get("order_note", "visible_in_input_svg")) for row in order] + additions
    combined.sort(key=lambda r: float(r["y"]))
    return combined


def beautify_tree_svg(svg_text, tip_rows, meta):
    out = svg_text
    out = re.sub(r'<g id="hoverHolder"[^>]*/>', "", out)
    out = re.sub(r'<g id="tooltipHolder".*?</g>', "", out, flags=re.S)
    out = re.sub(r'stroke="#000000"', 'stroke="#6B7280"', out)
    out = re.sub(r'stroke="#999999"', 'stroke="#9CA3AF"', out)
    out = re.sub(r'font-size="20"', 'font-size="22"', out)
    out = re.sub(r'font-family="Arial"', 'font-family="Arial, Helvetica, sans-serif"', out)

    by_label = {r["label"]: r for r in tip_rows}
    for label in sorted(by_label, key=len, reverse=True):
        habitat = habitat_for(label, meta)
        color = KEY.get(label, (None, None))[1] or HABITAT_COLORS.get(habitat, "#111827")
        weight = "700" if label in KEY or habitat in {"IS", "AS", "ES", "NS", "Unknown_ref"} else "400"
        pattern = (
            r'(<text\b(?=[^>]*\bx="1005")(?=[^>]*\by="' + re.escape(str(int(by_label[label]["y"]))) +
            r'")[^>]*?)fill="#000000"([^>]*>)' + re.escape(label) + r"(</text>)"
        )
        repl = r'\1fill="' + color + r'"\2' + esc(label) + r'\3'
        out = re.sub(pattern, repl, out)
        if weight == "700":
            pattern2 = (
                r'(<text\b(?=[^>]*\bx="1005")(?=[^>]*\by="' + re.escape(str(int(by_label[label]["y"]))) +
                r'")[^>]*?)font-weight=""([^>]*>)'
            )
            out = re.sub(pattern2, r'\1font-weight="700"\2', out)

    legend = [
        '<g id="codex-clean-legend" font-family="Arial, Helvetica, sans-serif">',
        '<rect x="1075" y="24" width="230" height="94" fill="#FFFFFF" fill-opacity="0.88" stroke="#D1D5DB" stroke-width="1"/>',
        '<text x="1090" y="45" font-size="18" font-weight="700" fill="#111827">Habitat / key</text>',
    ]
    lx, ly = 1090, 65
    for i, (label, color) in enumerate([
        ("IS", HABITAT_COLORS["IS"]),
        ("AS", HABITAT_COLORS["AS"]),
        ("ES", HABITAT_COLORS["ES"]),
        ("NS", HABITAT_COLORS["NS"]),
        ("Unknown ref", HABITAT_COLORS["Unknown_ref"]),
        ("Reference", HABITAT_COLORS["Reference"]),
    ]):
        x = lx + (i % 2) * 104
        y = ly + (i // 2) * 18
        legend.append(f'<rect x="{x}" y="{y-9}" width="12" height="10" fill="{color}"/>')
        legend.append(f'<text x="{x+18}" y="{y}" font-size="12" fill="#374151">{esc(label)}</text>')
    legend.append("</g>")
    out = out.replace("</svg>", "\n".join(legend) + "\n</svg>")
    return out


def arrow_points(x1, x2, y, h, strand):
    if x2 < x1:
        x1, x2 = x2, x1
    width = max(x2 - x1, 5)
    head = min(10, width * 0.42)
    y1, y2, ym = y - h / 2, y + h / 2, y
    if strand == "-":
        pts = [(x1, ym), (x1 + head, y1), (x2, y1), (x2, y2), (x1 + head, y2)]
    else:
        pts = [(x1, y1), (x2 - head, y1), (x2, ym), (x2 - head, y2), (x1, y2)]
    return " ".join(f"{a:.1f},{b:.1f}" for a, b in pts)


def make_arrow_svg(order, meta, manifest_rows):
    rows_by_label = {}
    for row in manifest_rows:
        rows_by_label.setdefault(row["tree_label"], []).append(row)
    augment_r2111_neighborhood(rows_by_label)
    mag_order = [r["label"] for r in order if r["label"] in rows_by_label]
    row_h = 18
    top = 120
    left = 32
    label_x = 76
    arrow_x = 610
    arrow_w = 800
    evidence_x = arrow_x + arrow_w + 28
    width = evidence_x + 330
    height = top + len(mag_order) * row_h + 90

    vals = []
    for label in mag_order:
        for r in rows_by_label[label]:
            try:
                vals.extend([int(float(r["relative_bp_start"])), int(float(r["relative_bp_end"]))])
            except Exception:
                pass
    rel_min = min(vals + [-12000])
    rel_max = max(vals + [14500])

    def sx(rel):
        return arrow_x + (rel - rel_min) / (rel_max - rel_min) * arrow_w

    chunks = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        '<text x="32" y="32" font-family="Arial, Helvetica, sans-serif" font-size="21" font-weight="700" fill="#111827">Small AioA/IdrA MAG gene-island arrows</text>',
        '<text x="32" y="56" font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#4B5563">Rows follow the supplied iTOL SVG top-to-bottom order. Only MAG rows with local gene-neighborhood/CD-search evidence are shown.</text>',
        '<text x="32" y="88" font-family="Arial, Helvetica, sans-serif" font-size="11" font-weight="700" fill="#111827">Sequence</text>',
        f'<text x="{arrow_x}" y="88" font-family="Arial, Helvetica, sans-serif" font-size="11" font-weight="700" fill="#111827">Gene island / neighborhood</text>',
        f'<text x="{evidence_x}" y="88" font-family="Arial, Helvetica, sans-serif" font-size="11" font-weight="700" fill="#111827">Prediction</text>',
        f'<line x1="{sx(0):.1f}" y1="{top-14}" x2="{sx(0):.1f}" y2="{height-70}" stroke="#94A3B8" stroke-width="0.8" stroke-dasharray="3 4"/>',
    ]
    lx = 32
    for label, color in [
        ("IS", HABITAT_COLORS["IS"]),
        ("AS", HABITAT_COLORS["AS"]),
        ("ES", HABITAT_COLORS["ES"]),
        ("NS", HABITAT_COLORS["NS"]),
        ("A", GENE_COLORS["AioA/IdrA-related A"]),
        ("B", GENE_COLORS["IdrB/AioB-related small subunit"]),
        ("P-like", GENE_COLORS["P-like"]),
        ("AioB", GENE_COLORS["canonical AioB"]),
        ("other", GENE_COLORS["other annotated CDS"]),
    ]:
        chunks.append(f'<rect x="{lx}" y="99" width="12" height="8" fill="{color}" stroke="#374151" stroke-width="0.3"/>')
        chunks.append(f'<text x="{lx+16}" y="107" font-family="Arial, Helvetica, sans-serif" font-size="10" fill="#374151">{esc(label)}</text>')
        lx += 76

    for i, label in enumerate(mag_order):
        y = top + i * row_h
        habitat = habitat_for(label, meta)
        hcolor = KEY.get(label, (None, None))[1] or HABITAT_COLORS.get(habitat, "#9CA3AF")
        rows = rows_by_label[label]
        display = meta.get(label, {}).get("display_label", label)
        if len(display) > 70:
            display = display[:67] + "..."
        chunks.append(f'<rect x="{left}" y="{y-5.5:.1f}" width="26" height="11" rx="1.5" fill="{hcolor}"/>')
        chunks.append(f'<text x="{label_x}" y="{y+3.3:.1f}" font-family="Arial, Helvetica, sans-serif" font-size="10" font-weight="{"700" if label in KEY else "400"}" fill="{hcolor}">{esc(display)}</text>')
        chunks.append(f'<line x1="{arrow_x}" y1="{y}" x2="{arrow_x+arrow_w}" y2="{y}" stroke="#E5E7EB" stroke-width="0.4"/>')
        for r in rows:
            try:
                x1 = sx(int(float(r["relative_bp_start"])))
                x2 = sx(int(float(r["relative_bp_end"])))
            except Exception:
                continue
            plot = r["plot_label"]
            fill = GENE_COLORS.get(plot, GENE_COLORS["other annotated CDS"])
            stroke = "#111827" if plot != "other annotated CDS" else "#94A3B8"
            chunks.append(f'<polygon points="{arrow_points(x1, x2, y, 10, r.get("strand", "+"))}" fill="{fill}" stroke="{stroke}" stroke-width="0.45"/>')
        pred = prediction_for_rows(rows)
        chunks.append(f'<text x="{evidence_x}" y="{y+3.3:.1f}" font-family="Arial, Helvetica, sans-serif" font-size="9.2" fill="{pred[1]}">{esc(pred[0])}</text>')
        if label in KEY:
            chunks.append(f'<text x="{arrow_x-136}" y="{y+3.3:.1f}" font-family="Arial, Helvetica, sans-serif" font-size="9" font-weight="700" fill="{hcolor}">{esc(KEY[label][0])}</text>')

    chunks.append(f'<text x="32" y="{height-32}" font-family="Arial, Helvetica, sans-serif" font-size="10" fill="#64748B">A = focal AioA/IdrA-related molybdopterin oxidoreductase; B = IdrB/AioB-related small subunit; P-like = two P-like CDSs identified by joint reference search.</text>')
    chunks.append("</svg>")
    return "\n".join(chunks), mag_order


def prediction_for_rows(rows):
    labels = [r.get("plot_label", "") for r in rows]
    has_aio_b = "canonical AioB" in labels
    has_b = "IdrB/AioB-related small subunit" in labels
    p_count = sum(1 for x in labels if x == "P-like")
    if has_b and p_count >= 2:
        return "IdrA-PP gene island: strong", "#146C43"
    if has_b and p_count == 1:
        return "IdrA-PP gene island: partial", "#B45309"
    if has_b:
        return "B-linked, P-like missing", "#B45309"
    if has_aio_b:
        return "canonical Aio-like", "#7C3AED"
    if "CD-search only A" in labels:
        return "focal A only, neighborhood missing", "#6B7280"
    return "unresolved/truncated", "#6B7280"


def augment_r2111_neighborhood(rows_by_label):
    label = "NS|R2111_N500_0-10_bin13-k141_4404301_3"
    rows = rows_by_label.get(label)
    if not rows:
        return
    labels = {row.get("plot_label", "") for row in rows}
    if "canonical AioB" in labels or "IdrB/AioB-related small subunit" in labels:
        return
    rows.insert(0, {
        "tree_label": label,
        "candidate_id": "R2111_N500_0-10_bin13-k141_4404301_3",
        "source": "remote reliable KOfam neighborhood",
        "gene_id": "R2111_N500_0-10_bin13-k141_4404301_2",
        "relative_gene_position": "-1",
        "relative_bp_start": "-1827",
        "relative_bp_end": "-1237",
        "strand": "+",
        "plot_label": "canonical AioB",
        "neighborhood_structure_class": "K08355-K08356 adjacent aio-like neighborhood",
        "strict_synteny_QC_pass": "not_available",
        "annotation": "K08355 arsenite oxidase small subunit; reliable KOfamScan (*)",
    })


def build_prediction_tables(order, meta, manifest_rows, strict_rows, hmm_rows, t640_rows):
    rows_by_label = {}
    for row in manifest_rows:
        rows_by_label.setdefault(row["tree_label"], []).append(row)
    augment_r2111_neighborhood(rows_by_label)
    strict_by_id = {r["candidate_id"]: r for r in strict_rows}
    hmm_by_id = {}
    for row in hmm_rows:
        pid = row["protein_id"]
        hmm_by_id[pid] = row
        hmm_by_id[strip_prefix(pid)] = row
    t640_by_norm = {r["normalized_id"]: r for r in t640_rows}

    pred_rows = []
    t640_small = []
    for rank, item in enumerate(order, start=1):
        label = item["label"]
        if label not in rows_by_label:
            continue
        cid = strip_prefix(label)
        rows = rows_by_label[label]
        pred, _ = prediction_for_rows(rows)
        strict = strict_by_id.get(cid, {})
        labels = [r.get("plot_label", "") for r in rows]
        pred_rows.append({
            "itol_order_rank": rank,
            "tree_label": label,
            "candidate_id": cid,
            "habitat": habitat_for(label, meta),
            "prediction": pred,
            "has_B_or_AioB_related": "yes" if ("IdrB/AioB-related small subunit" in labels or "canonical AioB" in labels) else "no",
            "P_like_arrow_count": str(sum(1 for x in labels if x == "P-like")),
            "strict_QC_pass": strict.get("QC_pass", "not_available"),
            "strict_neighborhood_class": strict.get("strict_neighborhood_class_revised") or strict.get("strict_neighborhood_class", "not_available"),
            "QC_fail_reasons": strict.get("QC_fail_reasons", ""),
        })
        hmm = hmm_by_id.get(cid, {})
        trow = t640_by_norm.get(cid, {})
        try:
            candidate_combined = float(hmm.get("combined_bitscore", "nan"))
        except ValueError:
            candidate_combined = float("nan")
        if trow:
            t640_status = "pass_T640_existing_allMAG_tblout"
            note = ""
        elif candidate_combined >= 640:
            t640_status = "candidate_HMM_pass_T640_not_in_allMAG_tblout"
            note = "Candidate HMM search passes score 640, but the all MAG genes.faa T640 table does not contain this normalized ID."
        else:
            t640_status = "not_in_existing_T640_tblout"
            note = "No all MAG T640 row was available for this small-tree candidate."
        t640_small.append({
            "itol_order_rank": rank,
            "tree_label": label,
            "candidate_id": cid,
            "habitat": habitat_for(label, meta),
            "combined_bitscore_48": hmm.get("combined_bitscore", ""),
            "idrA_bitscore": hmm.get("idrA_bitscore", ""),
            "aioA_bitscore": hmm.get("aioA_bitscore", ""),
            "idrA_minus_aioA": hmm.get("idrA_minus_aioA", ""),
            "T640_status": t640_status,
            "T640_combined_bitscore": trow.get("combined_bitscore", ""),
            "note": note,
        })
    return pred_rows, t640_small


def main():
    parser = argparse.ArgumentParser(description="Build small-tree iTOL-order gene island outputs.")
    parser.add_argument("--input-svg", type=Path, default=DEFAULT_INPUT_SVG)
    parser.add_argument("--treefile", type=Path, default=DEFAULT_TREEFILE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_META)
    parser.add_argument("--arrow-manifest", type=Path, default=DEFAULT_ARROW_MANIFEST)
    parser.add_argument("--hmm-scores", type=Path, default=DEFAULT_HMM_SCORES)
    parser.add_argument("--t640", type=Path, default=DEFAULT_T640)
    parser.add_argument("--strict-qc", type=Path, default=DEFAULT_STRICT_QC)
    args = parser.parse_args()

    args.out_dir.mkdir(exist_ok=True)
    svg_text = args.input_svg.read_text(errors="ignore")
    order = parse_svg_tip_order(svg_text)
    meta = {r["label"]: r for r in read_tsv(args.metadata)}
    manifest = read_tsv(args.arrow_manifest)
    manifest_labels = list(dict.fromkeys(row["tree_label"] for row in manifest))
    order = augment_order_from_treefile(order, manifest_labels, args.treefile)
    strict = read_tsv(args.strict_qc)
    hmm = read_tsv(args.hmm_scores)
    t640 = read_tsv(args.t640)

    order_path = args.out_dir / "itol_svg_tip_order.tsv"
    write_tsv(order_path, order, ["label", "x", "y", "order_note"])

    pretty_tree = beautify_tree_svg(svg_text, order, meta)
    pretty_path = args.out_dir / "c0J-IDrhpB2odQ8jRmUo-Q_beautified.svg"
    pretty_path.write_text(pretty_tree)

    arrow_svg, mag_order = make_arrow_svg(order, meta, manifest)
    arrow_path = args.out_dir / "small_tree_gene_island_arrows_iTOL_svg_order.svg"
    arrow_path.write_text(arrow_svg)

    pred_rows, t640_rows = build_prediction_tables(order, meta, manifest, strict, hmm, t640)
    pred_path = args.out_dir / "small_tree_idrA_PP_gene_island_predictions.tsv"
    write_tsv(pred_path, pred_rows, [
        "itol_order_rank", "tree_label", "candidate_id", "habitat", "prediction",
        "has_B_or_AioB_related", "P_like_arrow_count", "strict_QC_pass",
        "strict_neighborhood_class", "QC_fail_reasons",
    ])
    t640_path = args.out_dir / "small_tree_T640_threshold_check.tsv"
    write_tsv(t640_path, t640_rows, [
        "itol_order_rank", "tree_label", "candidate_id", "habitat",
        "combined_bitscore_48", "idrA_bitscore", "aioA_bitscore", "idrA_minus_aioA",
        "T640_status", "T640_combined_bitscore", "note",
    ])

    readme = args.out_dir / "README.txt"
    readme.write_text(
        "Small tree iTOL-order gene island outputs\n"
        "========================================\n"
        f"Input iTOL SVG: {args.input_svg}\n"
        f"Local small-tree folder: {BASE}\n"
        f"Treefile used to insert any missing MAG rows: {args.treefile}\n"
        f"Tip order parsed from SVG text y-coordinates: {order_path}\n\n"
        "Outputs:\n"
        f"- Beautified tree SVG: {pretty_path}\n"
        f"- Gene-island arrow SVG: {arrow_path}\n"
        f"- IdrA-PP/gene-island predictions: {pred_path}\n"
        f"- T640 threshold check for the small-tree candidates: {t640_path}\n\n"
        "Run note: /home/ps/... was not present in this local Codex environment, and hmmsearch/HMM model files were not available here. "
        "The T640 table therefore reuses the existing local tblout-derived HMM summaries and flags any small-tree candidate missing from those summaries.\n"
    )

    print(f"tips_in_svg_order={len(order)}")
    print(f"mag_rows_with_arrows={len(mag_order)}")
    print(f"predictions={len(pred_rows)}")
    print(f"t640_rows={len(t640_rows)}")
    print(pretty_path)
    print(arrow_path)
    print(pred_path)
    print(t640_path)


if __name__ == "__main__":
    main()
