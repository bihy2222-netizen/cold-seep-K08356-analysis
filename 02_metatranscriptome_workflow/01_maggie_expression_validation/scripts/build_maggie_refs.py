#!/usr/bin/env python3
import argparse
import csv
import hashlib
import os
from pathlib import Path


GENETIC_CODE = {
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


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_fasta(path):
    records = {}
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header.split()[0]] = (header, "".join(parts))
                header = line[1:]
                parts = []
            else:
                parts.append(line.strip())
    if header is not None:
        records[header.split()[0]] = (header, "".join(parts))
    return records


def wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def sha256_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def translate(seq):
    seq = seq.upper().replace("U", "T")
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa.append(GENETIC_CODE.get(codon, "X"))
    return "".join(aa)


def classify_row(row):
    final_class = row.get("Final_class", "")
    neigh = row.get("neighborhood_class", "")
    strict = row.get("strict_synteny_QC_pass", "").lower() == "yes"
    if strict:
        return "primary_strict_DIRM_like_IdrA"
    if final_class == "Canonical AioA clade":
        return "canonical_AioA"
    if "partial DIRM-like" in neigh:
        return "partial_IdrA_associated"
    if "Unknown DMSOR" in final_class:
        return "unknown_DMSOR"
    return "other_DMSOR_family"


def extract_record(row, prodigal_dir, faa_dir, prefix, serial):
    mag = row["MAG_ID"]
    protein_id = row["protein_ID"]
    cds_path = prodigal_dir / f"{mag}.cds.fna"
    faa_path = faa_dir / f"{mag}.cds.faa"
    result = {
        "ref_id": f"{prefix}{serial:03d}",
        "source_protein_ID": protein_id,
        "MAG_ID": mag,
        "cds_path": str(cds_path),
        "faa_path": str(faa_path),
        "extract_status": "PASS",
        "failure_reason": "",
    }
    if not cds_path.exists():
        result["extract_status"] = "FAIL"
        result["failure_reason"] = "missing_cds_file"
        return result, None
    if not faa_path.exists():
        result["extract_status"] = "FAIL"
        result["failure_reason"] = "missing_faa_file"
        return result, None
    cds_records = parse_fasta(cds_path)
    faa_records = parse_fasta(faa_path)
    if protein_id not in cds_records:
        result["extract_status"] = "FAIL"
        result["failure_reason"] = "protein_id_not_in_cds_fasta"
        return result, None
    if protein_id not in faa_records:
        result["extract_status"] = "FAIL"
        result["failure_reason"] = "protein_id_not_in_faa_fasta"
        return result, None
    cds_header, cds_seq = cds_records[protein_id]
    faa_header, faa_seq = faa_records[protein_id]
    cds_seq = cds_seq.upper()
    faa_seq = faa_seq.rstrip("*")
    tr = translate(cds_seq).rstrip("*")
    if tr == faa_seq:
        translation_status = "PASS_exact"
    elif tr.startswith(faa_seq) or faa_seq.startswith(tr):
        translation_status = "PASS_prefix_or_truncation_match"
    else:
        translation_status = "FAIL_translation_mismatch"
    result.update({
        "source_cds_header": cds_header,
        "source_faa_header": faa_header,
        "cds_length_nt": str(len(cds_seq)),
        "protein_length_aa": str(len(faa_seq)),
        "translated_length_aa": str(len(tr)),
        "translation_status": translation_status,
        "cds_sha256": sha256_text(cds_seq),
        "protein_sha256": sha256_text(faa_seq),
    })
    fasta_header = (
        f"{result['ref_id']} source={protein_id} MAG={mag} "
        f"class={classify_row(row)} final_class={row.get('Final_class','NA')} "
        f"neighborhood={row.get('neighborhood_class','NA').replace(' ', '_')}"
    )
    record = {
        "ref_id": result["ref_id"],
        "cds_header": fasta_header,
        "faa_header": fasta_header,
        "cds_seq": cds_seq,
        "faa_seq": faa_seq,
        "translation_status": translation_status,
    }
    return result, record


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_fasta(path, records, seq_key):
    with open(path, "w") as fh:
        for rec in records:
            fh.write(f">{rec[seq_key.replace('_seq', '_header')]}\n{wrap(rec[seq_key])}\n")


def duplicate_audit(records):
    rows = []
    seen_cds = {}
    seen_faa = {}
    for rec in records:
        cds_hash = sha256_text(rec["cds_seq"])
        faa_hash = sha256_text(rec["faa_seq"])
        rows.append({
            "ref_id": rec["ref_id"],
            "cds_sha256": cds_hash,
            "protein_sha256": faa_hash,
            "duplicate_cds_of": seen_cds.get(cds_hash, ""),
            "duplicate_protein_of": seen_faa.get(faa_hash, ""),
            "exact_duplicate_status": "DUPLICATE" if cds_hash in seen_cds or faa_hash in seen_faa else "UNIQUE",
        })
        seen_cds.setdefault(cds_hash, rec["ref_id"])
        seen_faa.setdefault(faa_hash, rec["ref_id"])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--neighborhood", required=True)
    ap.add_argument("--prodigal-dir", required=True)
    ap.add_argument("--faa-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(args.evidence)
    prodigal_dir = Path(args.prodigal_dir)
    faa_dir = Path(args.faa_dir)

    strict_rows = [r for r in rows if r.get("strict_synteny_QC_pass", "").lower() == "yes"]
    comp_rows = [
        r for r in rows
        if r.get("strict_synteny_QC_pass", "").lower() == "yes"
        or r.get("Final_class", "") in {"Canonical AioA clade", "Unknown DMSOR clade"}
        or "partial DIRM-like" in r.get("neighborhood_class", "")
    ]

    audit_rows = []
    strict_records = []
    for i, row in enumerate(strict_rows, 1):
        audit, rec = extract_record(row, prodigal_dir, faa_dir, "IdrA_HC", i)
        audit.update({
            "reference_set": "high_confidence_IdrA",
            "Final_class": row.get("Final_class", ""),
            "Confidence": row.get("Confidence", ""),
            "neighborhood_class": row.get("neighborhood_class", ""),
            "B_related_gene": row.get("B_related_gene", ""),
            "P_like_count": row.get("P_like_count", ""),
            "P_like_genes": row.get("P_like_genes", ""),
            "IdrA_bitscore": row.get("IdrA_bitscore", ""),
            "AioA_bitscore": row.get("AioA_bitscore", ""),
            "combined_bitscore": row.get("combined_bitscore", ""),
        })
        audit_rows.append(audit)
        if rec:
            strict_records.append(rec)

    comp_audit = []
    comp_records = []
    counters = {}
    for row in comp_rows:
        cls = classify_row(row)
        counters[cls] = counters.get(cls, 0) + 1
        prefix = {
            "primary_strict_DIRM_like_IdrA": "IdrA_HC",
            "partial_IdrA_associated": "IdrA_PART",
            "canonical_AioA": "AioA_CAN",
            "unknown_DMSOR": "DMSOR_UNK",
        }.get(cls, "DMSOR_OTHER")
        audit, rec = extract_record(row, prodigal_dir, faa_dir, prefix, counters[cls])
        audit.update({
            "reference_set": "competitive_DMSOR",
            "reference_class": cls,
            "Final_class": row.get("Final_class", ""),
            "Confidence": row.get("Confidence", ""),
            "neighborhood_class": row.get("neighborhood_class", ""),
            "IdrA_bitscore": row.get("IdrA_bitscore", ""),
            "AioA_bitscore": row.get("AioA_bitscore", ""),
            "combined_bitscore": row.get("combined_bitscore", ""),
        })
        comp_audit.append(audit)
        if rec:
            comp_records.append(rec)

    write_fasta(outdir / "high_confidence_IdrA_CDS.fna", strict_records, "cds_seq")
    write_fasta(outdir / "high_confidence_IdrA_proteins.faa", strict_records, "faa_seq")
    write_fasta(outdir / "competitive_DMSOR_CDS.fna", comp_records, "cds_seq")
    write_fasta(outdir / "competitive_DMSOR_proteins.faa", comp_records, "faa_seq")

    audit_fields = [
        "reference_set", "ref_id", "source_protein_ID", "MAG_ID", "Final_class", "Confidence",
        "neighborhood_class", "B_related_gene", "P_like_count", "P_like_genes",
        "IdrA_bitscore", "AioA_bitscore", "combined_bitscore",
        "cds_length_nt", "protein_length_aa", "translated_length_aa",
        "translation_status", "extract_status", "failure_reason",
        "cds_sha256", "protein_sha256", "cds_path", "faa_path", "source_cds_header", "source_faa_header",
    ]
    write_tsv(outdir / "02_high_confidence_IdrA_locus_audit.tsv", audit_rows, audit_fields)
    write_tsv(outdir / "02_high_confidence_IdrA_translation_check.tsv", audit_rows, audit_fields)
    write_tsv(outdir / "02_reference_manifest.tsv", comp_audit, [
        "reference_set", "reference_class", "ref_id", "source_protein_ID", "MAG_ID", "Final_class",
        "Confidence", "neighborhood_class", "IdrA_bitscore", "AioA_bitscore", "combined_bitscore",
        "cds_length_nt", "protein_length_aa", "translation_status", "extract_status",
        "failure_reason", "cds_sha256", "protein_sha256", "cds_path", "faa_path",
    ])
    dup_rows = duplicate_audit(comp_records)
    write_tsv(outdir / "02_reference_duplicate_audit.tsv", dup_rows, [
        "ref_id", "cds_sha256", "protein_sha256", "duplicate_cds_of", "duplicate_protein_of", "exact_duplicate_status",
    ])
    missing_classes = ["NapA", "ArrA", "DmsA"]
    with open(outdir / "BLOCKED_REFERENCE_LIBRARY.md", "w") as fh:
        fh.write("# Maggie reference-library audit\n\n")
        fh.write(f"- High-confidence IdrA strict records extracted: {len(strict_records)} / {len(strict_rows)}\n")
        fh.write(f"- Competitive MAG-derived DMSOR records extracted: {len(comp_records)} / {len(comp_rows)}\n")
        fh.write("- Missing required external nucleotide CDS classes for accepted competitive mapping: "
                 + ", ".join(missing_classes) + "\n")
        fh.write("- Decision: do not start accepted Bowtie2 expression validation until reliable nucleotide CDS references for these classes are added. Protein-only references or HMMs are not valid Bowtie2 targets.\n")
    with open(outdir / "commands_used.sh", "a") as fh:
        fh.write("python3 scripts/build_maggie_refs.py --evidence \"$T\" --neighborhood \"$N\" --prodigal-dir \"$PROD\" --faa-dir \"$FAA\" --outdir \"$ROOT\"\n")

    print(f"strict_extracted={len(strict_records)}/{len(strict_rows)}")
    print(f"competitive_extracted={len(comp_records)}/{len(comp_rows)}")
    print("blocked_required_external_cds=NapA,ArrA,DmsA")


if __name__ == "__main__":
    main()
