#!/usr/bin/env python3
import csv
import os
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path


ROOT = Path(os.environ.get("ROOT", "/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra"))
SUMMARY = ROOT / "11_summary"
EXPR = ROOT / "10_expression"
BOWTIE = ROOT / "06_bowtie2"
REF_DIR = ROOT / "05_idra_reference"
FASTP_SUMMARY = SUMMARY / "fastp_summary.tsv"


def read_tsv(path):
    if not path.exists():
        return []
    with open(path, newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_lengths(path):
    lengths = {}
    name = None
    n = 0
    if not path.exists():
        return lengths
    with open(path, errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line)
        if name is not None:
            lengths[name] = n
    return lengths


def count_clean_fragments(row):
    r1 = Path(row.get("clean_R1", ""))
    if not r1.exists():
        return 0
    try:
        out = subprocess.check_output(["seqkit", "stats", "-T", str(r1)], text=True)
        lines = out.strip().splitlines()
        if len(lines) >= 2:
            header = lines[0].split("\t")
            vals = lines[1].split("\t")
            data = dict(zip(header, vals))
            return int(data.get("num_seqs", "0").replace(",", ""))
    except Exception:
        pass
    try:
        out = subprocess.check_output(["bash", "-lc", f"gzip -cd {r1} | awk 'END{{print NR/4}}'"], text=True)
        return int(float(out.strip()))
    except Exception:
        return 0


def idxstats_counts(path):
    counts = {}
    if not path.exists():
        return counts
    with open(path) as handle:
        for line in handle:
            ref, length, mapped, unmapped = line.rstrip("\n").split("\t")[:4]
            if ref != "*":
                counts[ref] = {"length": int(length), "mapped_alignments": int(mapped)}
    return counts


def depth_stats(path, lengths):
    values = defaultdict(list)
    if path.exists():
        with open(path) as handle:
            for line in handle:
                ref, pos, depth = line.rstrip("\n").split("\t")[:3]
                values[ref].append(int(depth))
    out = {}
    for ref, length in lengths.items():
        depths = values.get(ref, [])
        if len(depths) < length:
            depths = depths + [0] * (length - len(depths))
        if not depths:
            depths = [0] * length
        thirds = max(length // 3, 1)
        five = depths[:thirds]
        middle = depths[thirds:2 * thirds]
        three = depths[2 * thirds:]
        out[ref] = {
            "mean_depth": sum(depths) / length if length else 0,
            "median_depth": statistics.median(depths) if depths else 0,
            "max_depth": max(depths) if depths else 0,
            "breadth_1x": sum(d >= 1 for d in depths) / length if length else 0,
            "breadth_3x": sum(d >= 3 for d in depths) / length if length else 0,
            "breadth_5x": sum(d >= 5 for d in depths) / length if length else 0,
            "breadth_10x": sum(d >= 10 for d in depths) / length if length else 0,
            "five_prime_breadth_1x": sum(d >= 1 for d in five) / len(five) if five else 0,
            "middle_breadth_1x": sum(d >= 1 for d in middle) / len(middle) if middle else 0,
            "three_prime_breadth_1x": sum(d >= 1 for d in three) / len(three) if three else 0,
            "zero_coverage_fraction": sum(d == 0 for d in depths) / length if length else 1,
        }
    return out


def evidence_level(mapq20_fragments, breadth_1x, regions_hit, competitive_note):
    if mapq20_fragments >= 5 and breadth_1x >= 0.20 and regions_hit >= 2 and competitive_note != "non_idrA_preferred":
        return "High confidence"
    if mapq20_fragments >= 3 and breadth_1x >= 0.10 and competitive_note != "non_idrA_preferred":
        return "Moderate confidence"
    if mapq20_fragments >= 1:
        return "Low confidence"
    return "Not detected"


def main():
    SUMMARY.mkdir(parents=True, exist_ok=True)
    fastp_rows = read_tsv(FASTP_SUMMARY)
    clean_fragments = {(r["sample"], r["run"]): count_clean_fragments(r) for r in fastp_rows}
    sample_by_run = {r["run"]: r["sample"] for r in fastp_rows}

    ref_fastas = {
        "strict": REF_DIR / "idra_strict_20_CDS.fna",
        "extended": REF_DIR / "idra_extended_CDS.fna",
        "competitive": REF_DIR / "DMSOR_competitive_CDS.fna",
    }
    rows = []
    for idxstats in sorted(EXPR.glob("*.idxstats.tsv")):
        parts = idxstats.name.split(".")
        if len(parts) < 3:
            continue
        run, panel = parts[0], parts[1]
        sample = sample_by_run.get(run, "")
        lengths = fasta_lengths(ref_fastas.get(panel, Path()))
        idx = idxstats_counts(idxstats)
        depth = depth_stats(EXPR / f"{run}.{panel}.MAPQ20.depth.tsv", lengths)
        clean = clean_fragments.get((sample, run), 0)
        for ref, length in lengths.items():
            mapped_alignments = idx.get(ref, {}).get("mapped_alignments", 0)
            # paired-end alignments count reads; approximate fragments as read alignments / 2
            raw_fragments = mapped_alignments / 2
            ds = depth.get(ref, {})
            mapq20_fragments = (ds.get("mean_depth", 0) * length) / 150 / 2 if length else 0
            fpm = (mapq20_fragments / clean * 1_000_000) if clean else 0
            rpkm = (mapq20_fragments * 1_000_000_000 / (length * clean)) if clean and length else 0
            regions_hit = sum([
                ds.get("five_prime_breadth_1x", 0) > 0,
                ds.get("middle_breadth_1x", 0) > 0,
                ds.get("three_prime_breadth_1x", 0) > 0,
            ])
            rows.append({
                "sample": sample,
                "run": run,
                "reference_panel": panel,
                "reference_gene": ref,
                "CDS_length_bp": length,
                "clean_fragments": clean,
                "raw_mapped_fragments_approx": f"{raw_fragments:.3f}",
                "MAPQ20_fragments_depth_approx": f"{mapq20_fragments:.3f}",
                "IdrA_FPM": f"{fpm:.6f}",
                "IdrA_RPKM": f"{rpkm:.6f}",
                "mean_depth": f"{ds.get('mean_depth', 0):.6f}",
                "median_depth": f"{ds.get('median_depth', 0):.6f}",
                "max_depth": f"{ds.get('max_depth', 0):.6f}",
                "breadth_1x": f"{ds.get('breadth_1x', 0):.6f}",
                "breadth_3x": f"{ds.get('breadth_3x', 0):.6f}",
                "breadth_5x": f"{ds.get('breadth_5x', 0):.6f}",
                "breadth_10x": f"{ds.get('breadth_10x', 0):.6f}",
                "five_prime_breadth_1x": f"{ds.get('five_prime_breadth_1x', 0):.6f}",
                "middle_breadth_1x": f"{ds.get('middle_breadth_1x', 0):.6f}",
                "three_prime_breadth_1x": f"{ds.get('three_prime_breadth_1x', 0):.6f}",
                "regions_with_coverage": regions_hit,
                "panel_mapping_note": "panel_FPM_RPKM_not_whole_transcriptome_TPM",
                "expression_evidence_level": evidence_level(mapq20_fragments, ds.get("breadth_1x", 0), regions_hit, ""),
            })

    fields = list(rows[0].keys()) if rows else [
        "sample", "run", "reference_panel", "reference_gene", "CDS_length_bp",
        "clean_fragments", "IdrA_FPM", "IdrA_RPKM", "mean_depth", "breadth_1x",
        "expression_evidence_level",
    ]
    out = SUMMARY / "idra_bowtie_expression.tsv"
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out} rows={len(rows)}")


if __name__ == "__main__":
    main()
