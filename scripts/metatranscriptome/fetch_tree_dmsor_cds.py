#!/usr/bin/env python3
import argparse
import csv
import hashlib
import time
import urllib.parse
import urllib.request
from pathlib import Path


CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def parse_fasta(path):
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    yield header, "".join(parts)
                header = line[1:]
                parts = []
            else:
                parts.append(line.strip())
        if header:
            yield header, "".join(parts)


def wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def translate(seq):
    seq = seq.upper().replace("U", "T")
    return "".join(CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def sha(seq):
    return hashlib.sha256(seq.encode()).hexdigest()


def infer_ref(header, seq, source):
    if header.startswith(("NapA|", "DmsA|", "ArxA|")):
        bits = header.split("|", 2)
        if len(bits) >= 2:
            cls, acc = bits[0], bits[1].split()[0]
            return cls, acc, header, seq, source
    if "__ArrA__" in header or "__ArxA__" in header:
        bits = header.split()
        parts = bits[0].split("__")
        if len(parts) >= 3:
            return parts[1], parts[2], header, seq, source
    return None


def fetch_cds(acc):
    query = urllib.parse.urlencode({
        "db": "protein",
        "id": acc,
        "rettype": "fasta_cds_na",
        "retmode": "text",
    })
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode()


def parse_single_fasta_text(text):
    lines = [x.strip() for x in text.splitlines() if x.strip()]
    if not lines or not lines[0].startswith(">"):
        return "", ""
    return lines[0][1:], "".join(lines[1:]).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--broad-faa", required=True)
    ap.add_argument("--arr-faa", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--per-class", type=int, default=6)
    args = ap.parse_args()

    out = Path(args.outdir)
    rows = []
    records = []
    candidates = []
    for path in [Path(args.broad_faa), Path(args.arr_faa)]:
        for header, seq in parse_fasta(path):
            rec = infer_ref(header, seq, str(path))
            if rec and rec[0] in {"NapA", "ArrA", "DmsA"}:
                candidates.append(rec)

    accepted = {"NapA": 0, "ArrA": 0, "DmsA": 0}
    seen_acc = set()
    for cls, acc, header, protein, source in candidates:
        if accepted[cls] >= args.per_class or acc in seen_acc:
            continue
        seen_acc.add(acc)
        status = "FAIL"
        reason = ""
        cds_header = ""
        cds_seq = ""
        try:
            text = fetch_cds(acc)
            cds_header, cds_seq = parse_single_fasta_text(text)
            if not cds_seq:
                reason = "no_fasta_cds_na_returned"
            else:
                translated = translate(cds_seq).rstrip("*")
                protein_clean = protein.rstrip("*")
                if translated == protein_clean:
                    status = "PASS_exact"
                elif translated in protein_clean or protein_clean in translated:
                    status = "PASS_contained_or_partial"
                else:
                    reason = "translation_mismatch"
        except Exception as exc:
            reason = f"fetch_error:{type(exc).__name__}:{exc}"
        row = {
            "class": cls,
            "accession": acc,
            "source_tree_header": header,
            "source_tree_file": source,
            "fetch_status": status,
            "failure_reason": reason,
            "cds_header": cds_header,
            "cds_length_nt": len(cds_seq),
            "protein_length_aa": len(protein.rstrip("*")),
            "cds_sha256": sha(cds_seq) if cds_seq else "",
            "protein_sha256": sha(protein.rstrip("*")),
        }
        rows.append(row)
        if status.startswith("PASS"):
            accepted[cls] += 1
            ref_id = f"EXT_{cls}_{accepted[cls]:02d}"
            records.append({
                "ref_id": ref_id,
                "class": cls,
                "accession": acc,
                "cds_header": f"{ref_id} class={cls} accession={acc} source=tree_reference",
                "faa_header": f"{ref_id} class={cls} accession={acc} source=tree_reference",
                "cds_seq": cds_seq,
                "faa_seq": protein.rstrip("*"),
                "fetch_row": row,
            })
        if all(v >= args.per_class for v in accepted.values()):
            break
        time.sleep(0.35)

    fields = ["class", "accession", "source_tree_header", "source_tree_file", "fetch_status",
              "failure_reason", "cds_header", "cds_length_nt", "protein_length_aa", "cds_sha256", "protein_sha256"]
    with open(out / "02_external_tree_DMSOR_CDS_fetch.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with open(out / "external_tree_DMSOR_CDS.fna", "w") as fna, open(out / "external_tree_DMSOR_proteins.faa", "w") as faa:
        for r in records:
            fna.write(f">{r['cds_header']}\n{wrap(r['cds_seq'])}\n")
            faa.write(f">{r['faa_header']}\n{wrap(r['faa_seq'])}\n")

    with open(out / "02_external_tree_DMSOR_reference_manifest.tsv", "w", newline="") as fh:
        fields2 = ["ref_id", "reference_class", "accession", "source_tree_header", "cds_length_nt",
                   "protein_length_aa", "translation_status", "cds_sha256", "protein_sha256"]
        writer = csv.DictWriter(fh, fieldnames=fields2, delimiter="\t")
        writer.writeheader()
        for r in records:
            fr = r["fetch_row"]
            writer.writerow({
                "ref_id": r["ref_id"],
                "reference_class": r["class"],
                "accession": r["accession"],
                "source_tree_header": fr["source_tree_header"],
                "cds_length_nt": fr["cds_length_nt"],
                "protein_length_aa": fr["protein_length_aa"],
                "translation_status": fr["fetch_status"],
                "cds_sha256": fr["cds_sha256"],
                "protein_sha256": fr["protein_sha256"],
            })
    print("accepted", accepted)
    if not all(v > 0 for v in accepted.values()):
        raise SystemExit("missing at least one required external class")


if __name__ == "__main__":
    main()
