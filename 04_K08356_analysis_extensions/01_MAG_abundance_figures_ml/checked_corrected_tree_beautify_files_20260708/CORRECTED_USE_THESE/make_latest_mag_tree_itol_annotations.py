from pathlib import Path


BASE = Path("/Users/catherine/Downloads/aioA 蛋白序列建树/MAG 补充 idra-tree")
TREE = BASE / "all_proteinstrimmed.fasta.treefile"
OUT = BASE / "latest_MAG_tree_beautify_files"


HABITAT_COLORS = {
    "IS": "#009E73",
    "AS": "#CC79A7",
    "ES": "#E69F00",
    "NS": "#0072B2",
}

REF_COLORS = {
    "MicrobiologySpectrum_AioA": "#E41A1C",
    "MicrobiologySpectrum_IdrA": "#6A3D9A",
    "ISMEJ_Fig5_AioA_or_AioA_like": "#00BFC4",
    "ISMEJ_Fig5_IdrA": "#4B0082",
    "ISMEJ_Fig5_unknown": "#F0E442",
    "ISMEJ_Fig5_other": "#7F7F7F",
    "Broad_DMSOR_reference": "#9E9E9E",
    "NarH_outgroup": "#111111",
}

TYPE_COLORS = {
    "Own_MAG_sample": "#D55E00",
    "Published_reference": "#6B7280",
    "Outgroup": "#111111",
}


def parse_newick_tips(text: str) -> list[str]:
    tips = []
    i = 0
    while i < len(text):
        if text[i] in "(,":
            j = i + 1
            if j < len(text) and text[j] == "(":
                i += 1
                continue
            k = j
            while k < len(text) and text[k] not in ":,);":
                k += 1
            label = text[j:k].strip()
            if label:
                tips.append(label)
            i = k
        else:
            i += 1
    return tips


def is_sample(label: str) -> bool:
    if "_bin" in label and "-k141_" in label:
        return True
    return label.startswith(("C1_", "C2_", "C3_", "R2111_", "S13_", "S14_", "S15_", "S1_", "SQ_", "SY"))


def habitat(label: str) -> str | None:
    if label.startswith(("SY365", "SY366", "SY368", "SY456", "SY457", "SY459")):
        return "IS"
    if label.startswith(("S13_", "S14_", "S15_", "S1_", "SQ_")):
        return "AS"
    if label.startswith(("C1_", "C2_", "C3_")):
        return "ES"
    if label.startswith("R2111_"):
        return "NS"
    return None


def reference_group(label: str) -> str | None:
    if label.endswith("_NarH") or "_NarH" in label:
        return "NarH_outgroup"
    if label.startswith("MicrobiologySpectrum_AioA|"):
        return "MicrobiologySpectrum_AioA"
    if label.startswith("MicrobiologySpectrum_IdrA|"):
        return "MicrobiologySpectrum_IdrA"
    if label.startswith("ISMEJFig5|"):
        if "|Arsenite_oxidase_clade|" in label:
            return "ISMEJ_Fig5_AioA_or_AioA_like"
        if "|Iodate_reductase_clade|" in label:
            return "ISMEJ_Fig5_IdrA"
        if "|Unknown_clade|" in label:
            return "ISMEJ_Fig5_unknown"
        return "ISMEJ_Fig5_other"
    if not is_sample(label):
        return "Broad_DMSOR_reference"
    return None


def write_colorstrip(path: Path, label: str, legend_title: str, colors: dict[str, str], rows: list[tuple[str, str]]):
    legend_keys = list(colors)
    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        f"DATASET_LABEL\t{label}",
        "COLOR\t#000000",
        f"LEGEND_TITLE\t{legend_title}",
        "LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_keys)),
        "LEGEND_COLORS\t" + "\t".join(colors[k] for k in legend_keys),
        "LEGEND_LABELS\t" + "\t".join(legend_keys),
        "DATA",
    ]
    for seq_id, group in rows:
        lines.append(f"{seq_id}\t{colors[group]}\t{group}")
    path.write_text("\n".join(lines) + "\n")


def write_tree_colors(path: Path, tips: list[str]):
    rows = []
    for label in tips:
        h = habitat(label)
        ref = reference_group(label)
        if h:
            color = HABITAT_COLORS[h]
            rows.append(f"{label}\tlabel\t{color}\tbold\t1")
            rows.append(f"{label}\tbranch\t{color}\tnormal\t3")
        elif ref:
            color = REF_COLORS[ref]
            width = "2.4" if ref != "Broad_DMSOR_reference" else "1.4"
            rows.append(f"{label}\tlabel\t{color}\tnormal\t1")
            rows.append(f"{label}\tbranch\t{color}\tnormal\t{width}")
    path.write_text("TREE_COLORS\nSEPARATOR TAB\nDATA\n" + "\n".join(rows) + "\n")


def main():
    OUT.mkdir(exist_ok=True)
    tips = parse_newick_tips(TREE.read_text())

    habitat_rows = [(tip, habitat(tip)) for tip in tips if habitat(tip)]
    reference_rows = [(tip, reference_group(tip)) for tip in tips if reference_group(tip)]
    type_rows = []
    for tip in tips:
        if habitat(tip):
            type_rows.append((tip, "Own_MAG_sample"))
        elif reference_group(tip) == "NarH_outgroup":
            type_rows.append((tip, "Outgroup"))
        else:
            type_rows.append((tip, "Published_reference"))

    write_colorstrip(
        OUT / "itol_COLORSTRIP_own_MAG_habitat.txt",
        "Own MAG habitat",
        "Own MAG habitat",
        HABITAT_COLORS,
        habitat_rows,
    )
    write_colorstrip(
        OUT / "itol_COLORSTRIP_reference_sources.txt",
        "Reference sources",
        "Reference sources",
        REF_COLORS,
        reference_rows,
    )
    write_colorstrip(
        OUT / "itol_COLORSTRIP_sequence_type.txt",
        "Sequence type",
        "Sequence type",
        TYPE_COLORS,
        type_rows,
    )
    write_tree_colors(OUT / "itol_TREE_COLORS_highlight_own_MAG_and_refs.txt", tips)
    (OUT / "latest_MAG_tree_for_iTOL.treefile").write_text(TREE.read_text())

    counts = {}
    for _, group in habitat_rows:
        counts[f"habitat:{group}"] = counts.get(f"habitat:{group}", 0) + 1
    for _, group in reference_rows:
        counts[f"reference:{group}"] = counts.get(f"reference:{group}", 0) + 1
    lines = ["category\tcount"] + [f"{k}\t{v}" for k, v in sorted(counts.items())]
    (OUT / "annotation_summary.tsv").write_text("\n".join(lines) + "\n")
    (OUT / "README_iTOL_upload_order.txt").write_text(
        "Upload latest_MAG_tree_for_iTOL.treefile to iTOL first, then drag these annotation files onto the tree:\n"
        "1. itol_COLORSTRIP_sequence_type.txt\n"
        "2. itol_COLORSTRIP_own_MAG_habitat.txt\n"
        "3. itol_COLORSTRIP_reference_sources.txt\n"
        "4. itol_TREE_COLORS_highlight_own_MAG_and_refs.txt\n\n"
        "Tip labels are unchanged. Own MAG/sample sequences are highlighted by habitat on labels and branches.\n"
    )
    print(OUT)


if __name__ == "__main__":
    main()
