from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


BASE = Path("/Users/catherine/Downloads/aioA 蛋白序列建树/MAG 补充 idra-tree")
WORK = BASE / "latest_MAG_tree_beautify_files"
TREE = BASE / "all_proteinstrimmed.fasta.treefile"
OUT_JSON = WORK / "MAG_tree_clade_statistics.json"
OUT_DETAIL = WORK / "MAG_tree_clade_detail.csv"
OUT_LISTS = WORK / "MAG_tree_clade_lists_like_tree11.tsv"


CLADE_ORDER = [
    "IdrA-associated",
    "canonical aioA-associated",
    "aioA-like-associated",
    "unknown AioA-like / uncertain DMSOR",
    "other DMSOR-family reference-associated",
]


class Node:
    def __init__(self) -> None:
        self.name = ""
        self.length = 0.0
        self.parent: Node | None = None
        self.children: list[Node] = []

    @property
    def is_leaf(self) -> bool:
        return not self.children


def parse_newick(text: str) -> Node:
    text = text.strip().rstrip(";")
    i = 0

    def parse_subtree() -> Node:
        nonlocal i
        node = Node()
        if text[i] == "(":
            i += 1
            while True:
                child = parse_subtree()
                child.parent = node
                node.children.append(child)
                if text[i] == ",":
                    i += 1
                    continue
                if text[i] == ")":
                    i += 1
                    break
        start = i
        while i < len(text) and text[i] not in ":,()":
            i += 1
        node.name = text[start:i].strip()
        if i < len(text) and text[i] == ":":
            i += 1
            start = i
            while i < len(text) and text[i] not in ",()":
                i += 1
            try:
                node.length = float(text[start:i])
            except ValueError:
                node.length = 0.0
        return node

    return parse_subtree()


def leaves(root: Node) -> list[Node]:
    stack = [root]
    out = []
    while stack:
        node = stack.pop()
        if node.is_leaf:
            out.append(node)
        else:
            stack.extend(reversed(node.children))
    return out


def is_mag(label: str) -> bool:
    return "_bin" in label and "-k141_" in label


def habitat(label: str) -> str:
    if label.startswith(("SY365", "SY366", "SY368", "SY456", "SY457", "SY459")):
        return "IS"
    if label.startswith(("S13_", "S14_", "S15_", "S1_", "SQ_")):
        return "AS"
    if label.startswith(("C1_", "C2_", "C3_")):
        return "ES"
    if label.startswith("R2111_"):
        return "NS"
    return "NA"


def sample_group(label: str) -> str:
    if label.startswith("R2111_"):
        return "_".join(label.split("_")[:2])
    if label.startswith("SQ_"):
        return "_".join(label.split("_")[:2])
    if label.startswith("S") and not label.startswith("SY"):
        return label.split("_", 1)[0]
    if label.startswith("C"):
        return label.split("_", 1)[0]
    if label.startswith("SY"):
        return label.split("-", 1)[0]
    return label.split("_", 1)[0]


def reference_type(label: str) -> str:
    if label.endswith("_NarH") or "_NarH" in label:
        return "NarH outgroup"
    if label.startswith("MicrobiologySpectrum_AioA|"):
        return "MicrobiologySpectrum AioA"
    if label.startswith("MicrobiologySpectrum_IdrA|"):
        return "MicrobiologySpectrum IdrA"
    if label.startswith("ISMEJFig5|Arsenite_oxidase_clade|"):
        return "ISMEJ Fig5 AioA"
    if label.startswith("ISMEJFig5|Iodate_reductase_clade|"):
        return "ISMEJ Fig5 IdrA"
    if label.startswith("ISMEJFig5|Unknown_clade|"):
        return "ISMEJ Fig5 unknown DMSOR"
    if "_AioA_Like" in label:
        return "AioA-like reference"
    if label.endswith("_AioA") or label.endswith("_aioA"):
        return "canonical AioA reference"
    if label.endswith("_Ahy"):
        return "Ahy/AioA-homolog reference"
    suffix = label.rsplit("_", 1)[-1]
    return f"{suffix} reference"


def clade_from_ref_type(ref_type: str) -> str:
    if "IdrA" in ref_type:
        return "IdrA-associated"
    if ref_type in {"MicrobiologySpectrum AioA", "ISMEJ Fig5 AioA", "canonical AioA reference"}:
        return "canonical aioA-associated"
    if ref_type == "AioA-like reference":
        return "aioA-like-associated"
    if ref_type in {"ISMEJ Fig5 unknown DMSOR", "Ahy/AioA-homolog reference", "Unka reference"}:
        return "unknown AioA-like / uncertain DMSOR"
    return "other DMSOR-family reference-associated"


def depth_map(root: Node) -> dict[Node, float]:
    depths = {root: 0.0}
    stack = [root]
    while stack:
        node = stack.pop()
        for child in node.children:
            depths[child] = depths[node] + child.length
            stack.append(child)
    return depths


def ancestors(node: Node) -> list[Node]:
    out = []
    while node is not None:
        out.append(node)
        node = node.parent
    return out


def distance(a: Node, b: Node, depths: dict[Node, float]) -> float:
    a_anc = {node: depths[node] for node in ancestors(a)}
    node = b
    while node not in a_anc:
        node = node.parent
        if node is None:
            return 0.0
    lca = node
    return depths[a] + depths[b] - 2 * depths[lca]


def main() -> None:
    root = parse_newick(TREE.read_text())
    tip_nodes = leaves(root)
    depths = depth_map(root)
    mag_nodes = [node for node in tip_nodes if is_mag(node.name)]
    ref_nodes = [node for node in tip_nodes if not is_mag(node.name)]

    rows = []
    for mag in mag_nodes:
        ranked = sorted(
            (
                (distance(mag, ref, depths), ref.name, reference_type(ref.name))
                for ref in ref_nodes
                if "NarH outgroup" not in reference_type(ref.name)
            ),
            key=lambda x: x[0],
        )
        nearest = ranked[0]
        nearest3 = ranked[:3]
        nearest5_clades = [clade_from_ref_type(item[2]) for item in ranked[:5]]
        clade_votes = Counter(nearest5_clades)
        assigned = clade_from_ref_type(nearest[2])
        rows.append(
            {
                "assigned_clade": assigned,
                "mag_id": mag.name,
                "habitat": habitat(mag.name),
                "sample_group": sample_group(mag.name),
                "nearest_reference": nearest[1],
                "nearest_reference_type": nearest[2],
                "tree_distance_to_nearest_reference": round(nearest[0], 6),
                "nearest_3_references": "; ".join(f"{r} [{t}]" for _, r, t in nearest3),
                "nearest_5_clade_votes": "; ".join(f"{k}:{v}" for k, v in clade_votes.most_common()),
                "display": f"{mag.name} | {habitat(mag.name)} | {sample_group(mag.name)}",
            }
        )

    rows.sort(key=lambda r: (CLADE_ORDER.index(r["assigned_clade"]), r["habitat"], r["sample_group"], r["mag_id"]))
    counts = Counter(row["assigned_clade"] for row in rows)
    habitat_counts = Counter((row["assigned_clade"], row["habitat"]) for row in rows)

    by_clade = defaultdict(list)
    for row in rows:
        by_clade[row["assigned_clade"]].append(row["display"])
    headers = [f"{clade} (n={counts[clade]})" for clade in CLADE_ORDER]
    max_len = max((len(by_clade[clade]) for clade in CLADE_ORDER), default=0)
    matrix = [headers]
    for idx in range(max_len):
        matrix.append([by_clade[clade][idx] if idx < len(by_clade[clade]) else "" for clade in CLADE_ORDER])

    with OUT_DETAIL.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    with OUT_LISTS.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows(matrix)

    OUT_JSON.write_text(
        json.dumps(
            {
                "matrix": matrix,
                "detail": rows,
                "counts": {clade: counts[clade] for clade in CLADE_ORDER},
                "habitat_counts": [
                    {"assigned_clade": clade, "habitat": hab, "n": n}
                    for (clade, hab), n in sorted(habitat_counts.items())
                ],
                "clade_order": CLADE_ORDER,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    print(f"MAG rows={len(rows)}")
    print(dict(counts))
    print(OUT_JSON)


if __name__ == "__main__":
    main()
