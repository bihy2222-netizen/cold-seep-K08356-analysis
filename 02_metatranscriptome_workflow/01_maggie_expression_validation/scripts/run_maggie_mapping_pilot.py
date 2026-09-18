#!/usr/bin/env python3
import argparse
import csv
import json
import math
import os
import shutil
import subprocess
from pathlib import Path


def run(cmd, log_path, cwd=None):
    with open(log_path, "a") as log:
        log.write("$ " + " ".join(map(str, cmd)) + "\n")
        log.flush()
        subprocess.run(cmd, cwd=cwd, check=True, stdout=log, stderr=subprocess.STDOUT)


def parse_fasta_lengths(path):
    lengths = {}
    header = None
    size = 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    lengths[header] = size
                header = line[1:].split()[0]
                size = 0
            else:
                size += len(line)
        if header:
            lengths[header] = size
    return lengths


def parse_idxstats(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            ref, length, mapped, unmapped = line.rstrip("\n").split("\t")
            if ref == "*":
                continue
            rows[ref] = {
                "length": int(length),
                "mapped_reads": int(mapped),
                "unmapped_reads_on_ref": int(unmapped),
            }
    return rows


def depth_stats(depth_path, lengths):
    covered = {k: 0 for k in lengths}
    total_depth = {k: 0 for k in lengths}
    first = {k: 0 for k in lengths}
    mid = {k: 0 for k in lengths}
    last = {k: 0 for k in lengths}
    with open(depth_path) as fh:
        for line in fh:
            ref, pos, dep = line.rstrip("\n").split("\t")[:3]
            if ref not in lengths:
                continue
            pos = int(pos)
            dep = int(dep)
            if dep > 0:
                covered[ref] += 1
            total_depth[ref] += dep
            length = lengths[ref]
            frac = pos / max(length, 1)
            if frac <= 1/3:
                first[ref] += dep
            elif frac <= 2/3:
                mid[ref] += dep
            else:
                last[ref] += dep
    stats = {}
    for ref, length in lengths.items():
        mean_depth = total_depth[ref] / length if length else 0.0
        breadth = covered[ref] / length if length else 0.0
        vals = [first[ref], mid[ref], last[ref]]
        nz = [v for v in vals if v > 0]
        uniformity = min(nz) / max(nz) if len(nz) == 3 and max(nz) else 0.0
        stats[ref] = {
            "covered_bases": covered[ref],
            "mean_depth": mean_depth,
            "coverage_breadth": breadth,
            "depth_5prime_third": first[ref],
            "depth_middle_third": mid[ref],
            "depth_3prime_third": last[ref],
            "coverage_uniformity_5_3": uniformity,
        }
    return stats


def parse_flagstat_json(path):
    try:
        data = json.loads(Path(path).read_text())
    except Exception:
        return {}
    qc = data.get("QC-passed reads", {})
    return {
        "total_reads": qc.get("total reads", 0),
        "mapped_reads_flagstat": qc.get("mapped reads", 0),
        "properly_paired_reads": qc.get("properly paired", 0),
    }


def count_bam_filters(bam, samtools, threads, log):
    def count(args):
        out = subprocess.check_output([samtools, "view", "-@", str(threads), "-c"] + args + [str(bam)])
        return int(out.decode().strip())
    return {
        "mapq20_reads": count(["-q", "20", "-F", "4"]),
        "mapq30_reads": count(["-q", "30", "-F", "4"]),
        "proper_pair_reads": count(["-f", "2", "-F", "4"]),
    }


def write_summary(out_path, sample, mode, idxstats, depths, bam_filters, flagstat, total_mapped_reads):
    total_counts = sum(v["mapped_reads"] for v in idxstats.values())
    denom_rpk = sum((v["mapped_reads"] / (v["length"] / 1000)) for v in idxstats.values() if v["length"]) or 0
    rows = []
    for ref, st in idxstats.items():
        length = st["length"]
        mapped = st["mapped_reads"]
        rpk = mapped / (length / 1000) if length else 0
        tpm = (rpk / denom_rpk * 1_000_000) if denom_rpk else 0
        cpm = (mapped / total_counts * 1_000_000) if total_counts else 0
        d = depths.get(ref, {})
        breadth = d.get("coverage_breadth", 0.0)
        mean_depth = d.get("mean_depth", 0.0)
        if mapped == 0:
            classification = "no transcriptional support"
        elif breadth >= 0.50 and mean_depth >= 1 and d.get("coverage_uniformity_5_3", 0) >= 0.20:
            is_idra = "IdrA_HC" in ref or "personalized_strict_IdrA" in ref
            classification = "robust IdrA-supported" if is_idra else "non-IdrA DMSOR-supported"
        elif breadth < 0.20:
            classification = "ambiguous DMSOR cross-mapping"
        else:
            classification = "weak_or_partial_transcriptional_support"
        row = {
            "sample": sample,
            "mode": mode,
            "reference_id": ref,
            "reference_length": length,
            "mapped_reads": mapped,
            "CPM": f"{cpm:.6f}",
            "TPM": f"{tpm:.6f}",
            "mean_depth": f"{mean_depth:.6f}",
            "coverage_breadth": f"{breadth:.6f}",
            "covered_bases": d.get("covered_bases", 0),
            "coverage_uniformity_5_3": f"{d.get('coverage_uniformity_5_3', 0.0):.6f}",
            "mapq20_reads_bam": bam_filters.get("mapq20_reads", 0),
            "mapq30_reads_bam": bam_filters.get("mapq30_reads", 0),
            "proper_pair_reads_bam": bam_filters.get("proper_pair_reads", 0),
            "total_mapped_reads_flagstat": flagstat.get("mapped_reads_flagstat", 0),
            "final_classification": classification,
        }
        rows.append(row)
    fields = list(rows[0].keys()) if rows else []
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--primary-ref", required=True)
    ap.add_argument("--competitive-ref", required=True)
    ap.add_argument("--primary-mode", default="primary_strict20")
    ap.add_argument("--competitive-mode", default="competitive_63")
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    root = Path(args.root)
    out = root / "maggie_mapping_pilot" / args.sample
    logs = out / "logs"
    refs = out / "refs"
    bam_dir = out / "bam"
    tables = out / "tables"
    igv = out / "IGV_files"
    for d in [logs, refs, bam_dir, tables, igv]:
        d.mkdir(parents=True, exist_ok=True)

    bowtie2 = shutil.which("bowtie2") or "bowtie2"
    bowtie2_build = shutil.which("bowtie2-build") or "bowtie2-build"
    samtools = shutil.which("samtools") or "samtools"
    log = logs / "maggie_mapping_pilot.log"

    modes = [
        (args.primary_mode, Path(args.primary_ref)),
        (args.competitive_mode, Path(args.competitive_ref)),
    ]
    all_rows = []
    for mode, ref in modes:
        ref_copy = refs / f"{mode}.fna"
        if not ref_copy.exists():
            shutil.copy2(ref, ref_copy)
        index_prefix = refs / mode
        if not (refs / f"{mode}.1.bt2").exists() and not (refs / f"{mode}.1.bt2l").exists():
            run([bowtie2_build, str(ref_copy), str(index_prefix)], log)
        bam = bam_dir / f"{args.sample}.{mode}.sorted.bam"
        if not bam.exists():
            cmd = (
                f"{bowtie2} -x {index_prefix} -1 {args.r1} -2 {args.r2} "
                f"--very-sensitive-local --no-mixed --no-discordant -k 10 -p {args.threads} "
                f"2> {logs / (mode + '.bowtie2.stderr.log')} | "
                f"{samtools} view -@ {args.threads} -bS - | "
                f"{samtools} sort -@ {args.threads} -o {bam} -"
            )
            run(["bash", "-lc", cmd], log)
        bai = Path(str(bam) + ".bai")
        if not bai.exists():
            run([samtools, "index", str(bam)], log)
        flag_json = tables / f"{args.sample}.{mode}.flagstat.json"
        run([samtools, "flagstat", "-O", "json", str(bam)], flag_json)
        idx = tables / f"{args.sample}.{mode}.idxstats.tsv"
        with open(idx, "w") as fh:
            subprocess.run([samtools, "idxstats", str(bam)], check=True, stdout=fh)
        depth = tables / f"{args.sample}.{mode}.depth.tsv"
        with open(depth, "w") as fh:
            subprocess.run([samtools, "depth", "-aa", str(bam)], check=True, stdout=fh)
        shutil.copy2(bam, igv / bam.name)
        shutil.copy2(bai, igv / bai.name)
        shutil.copy2(ref_copy, igv / ref_copy.name)
        lengths = parse_fasta_lengths(ref_copy)
        idxstats = parse_idxstats(idx)
        depths = depth_stats(depth, lengths)
        flagstat = parse_flagstat_json(flag_json)
        filters = count_bam_filters(bam, samtools, args.threads, log)
        rows = write_summary(
            tables / f"{args.sample}.{mode}.per_cds_summary.tsv",
            args.sample, mode, idxstats, depths, filters, flagstat, sum(x["mapped_reads"] for x in idxstats.values())
        )
        all_rows.extend(rows)

    combined = tables / f"{args.sample}.primary_vs_competitive_per_cds_summary.tsv"
    with open(combined, "w", newline="") as fh:
        fields = list(all_rows[0].keys()) if all_rows else []
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)

    report = out / f"{args.sample}.maggie_mapping_pilot_report.md"
    with open(report, "w") as fh:
        fh.write(f"# Maggie mapping pilot: {args.sample}\n\n")
        primary_count = len(parse_fasta_lengths(args.primary_ref))
        competitive_count = len(parse_fasta_lengths(args.competitive_ref))
        fh.write(f"- Primary reference: {args.primary_mode}, {primary_count} records.\n")
        fh.write(f"- Competitive reference: {args.competitive_mode}, {competitive_count} records.\n")
        fh.write("- Mapping: Bowtie2 very-sensitive-local, paired only, no mixed/no discordant, k=10.\n")
        fh.write("- Tables are in `tables/`; IGV BAM/BAI/reference files are in `IGV_files/`.\n")


if __name__ == "__main__":
    main()
