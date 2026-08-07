#!/usr/bin/env python3
import argparse
import csv
import math
import os
import re
from collections import defaultdict
from statistics import mean, median


HABITAT_ORDER = ["IS", "AS", "ES", "NS"]


def infer_habitat(sample):
    if sample.startswith("SY"):
        return "IS"
    if sample.startswith("SQ_"):
        return "AS"
    if sample.startswith(("S13", "S14", "S15")):
        return "ES"
    if sample.startswith("S3"):
        return "NS"
    if sample.startswith(("S1", "S2", "S4")):
        return "AS"
    if sample.startswith(("C", "ES_")):
        return "ES"
    if sample.startswith(("NS_", "R2111_")):
        return "NS"
    return "Unknown"


def read_matrix(path):
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        rows = []
        for row in reader:
            if not row:
                continue
            rows.append((row[0], [float(x) if x else 0.0 for x in row[1:]]))
    return samples, rows


def write_sample_map(samples, path):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Sample", "Habitat", "Rule"])
        for sample in samples:
            habitat = infer_habitat(sample)
            writer.writerow([sample, habitat, "inferred_from_sample_name"])


def read_sample_map(path, samples):
    mapping = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mapping[row["Sample"]] = row["Habitat"]
    missing = [s for s in samples if s not in mapping]
    if missing:
        raise SystemExit(f"Sample map missing {len(missing)} samples: {', '.join(missing[:8])}")
    return mapping


def rankdata(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg_rank
        i = j
    return ranks


def normal_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def chi2_sf(x, df):
    # Exact survival functions for df 1, 2, 3 are enough for these four-habitat tests.
    if df == 1:
        return math.erfc(math.sqrt(x / 2.0))
    if df == 2:
        return math.exp(-x / 2.0)
    if df == 3:
        z = math.sqrt(x)
        return math.erfc(z / math.sqrt(2.0)) + math.sqrt(2.0 / math.pi) * z * math.exp(-x / 2.0)
    # Wilson-Hilferty normal approximation fallback.
    z = ((x / df) ** (1.0 / 3.0) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
    return 1.0 - normal_cdf(z)


def kruskal(groups):
    data = []
    for group, vals in groups.items():
        for val in vals:
            data.append((group, val))
    n = len(data)
    if n == 0 or len(groups) < 2:
        return float("nan"), float("nan")
    ranks = rankdata([v for _, v in data])
    rank_sums = defaultdict(float)
    counts = defaultdict(int)
    for (group, _), rank in zip(data, ranks):
        rank_sums[group] += rank
        counts[group] += 1
    h = 12.0 / (n * (n + 1)) * sum(rank_sums[g] ** 2 / counts[g] for g in counts) - 3 * (n + 1)
    tie_counts = defaultdict(int)
    for _, val in data:
        tie_counts[val] += 1
    tie_term = sum(t ** 3 - t for t in tie_counts.values())
    denom = n ** 3 - n
    if denom and tie_term:
        correction = 1.0 - tie_term / denom
        if correction > 0:
            h /= correction
    return h, chi2_sf(max(h, 0.0), len(counts) - 1)


def mann_whitney_u(x, y):
    vals = [(v, 0) for v in x] + [(v, 1) for v in y]
    ranks = rankdata([v for v, _ in vals])
    r1 = sum(r for r, (_, group) in zip(ranks, vals) if group == 0)
    n1, n2 = len(x), len(y)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    mean_u = n1 * n2 / 2.0
    tie_counts = defaultdict(int)
    for v, _ in vals:
        tie_counts[v] += 1
    n = n1 + n2
    tie_term = sum(t ** 3 - t for t in tie_counts.values())
    var_u = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1)) if n > 1 else 0.0)
    if var_u <= 0:
        return u1, 1.0
    z = (abs(u1 - mean_u) - 0.5) / math.sqrt(var_u)
    p = 2.0 * (1.0 - normal_cdf(z))
    return u1, max(0.0, min(1.0, p))


def bh_adjust(pvals):
    n = len(pvals)
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    adj = [1.0] * n
    running = 1.0
    for rank, (idx, p) in reversed(list(enumerate(indexed, start=1))):
        running = min(running, p * n / rank)
        adj[idx] = min(1.0, running)
    return adj


def summarize(groups):
    out = {}
    for habitat in HABITAT_ORDER:
        vals = groups.get(habitat, [])
        out[habitat] = {
            "n": len(vals),
            "mean": mean(vals) if vals else float("nan"),
            "median": median(vals) if vals else float("nan"),
            "nonzero_n": sum(v > 0 for v in vals),
        }
    return out


def write_kw(path, results):
    with open(path, "w", newline="") as f:
        fieldnames = ["Feature", "Kruskal_H", "Kruskal_p"] + [
            f"{h}_{stat}" for h in HABITAT_ORDER for stat in ["n", "mean", "median", "nonzero_n"]
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for feature, groups in results.items():
            h, p = kruskal(groups)
            stats = summarize(groups)
            row = {"Feature": feature, "Kruskal_H": f"{h:.8g}", "Kruskal_p": f"{p:.8g}"}
            for habitat in HABITAT_ORDER:
                for stat, val in stats[habitat].items():
                    row[f"{habitat}_{stat}"] = f"{val:.8g}" if isinstance(val, float) else val
            writer.writerow(row)


def write_pairwise(path, results):
    rows = []
    for feature, groups in results.items():
        habitats = [h for h in HABITAT_ORDER if h in groups]
        raw = []
        start = len(rows)
        for i, h1 in enumerate(habitats):
            for h2 in habitats[i + 1:]:
                u, p = mann_whitney_u(groups[h1], groups[h2])
                rows.append({
                    "Feature": feature, "Group1": h1, "Group2": h2,
                    "U": f"{u:.8g}", "p": p,
                })
                raw.append(p)
        adj = bh_adjust(raw)
        for offset, p_adj in enumerate(adj):
            rows[start + offset]["p_adj_BH"] = p_adj
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Feature", "Group1", "Group2", "U", "p", "p_adj_BH"], delimiter="\t")
        writer.writeheader()
        for row in rows:
            row = dict(row)
            row["p"] = f"{row['p']:.8g}"
            row["p_adj_BH"] = f"{row['p_adj_BH']:.8g}"
            writer.writerow(row)


def clean_mag_id(label):
    label = re.sub(r"\s+", "_", label.strip())
    label = label.replace("_bin", "_bin")
    m = re.search(r"(.+?_bin\d+)", label)
    return m.group(1) if m else None


def read_feature_map(path):
    feature_to_mags = defaultdict(list)
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "MAG" not in reader.fieldnames or "Feature" not in reader.fieldnames:
            raise SystemExit(f"{path} must contain columns: MAG and Feature")
        for row in reader:
            mag = row["MAG"].strip()
            feature = row["Feature"].strip()
            if mag and feature:
                feature_to_mags[feature].append(mag)
    return feature_to_mags


def write_long(path, samples, sample_to_habitat, feature_values):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Feature", "Sample", "Habitat", "Abundance", "log10_abundance_plus_1"])
        for feature, values in feature_values.items():
            for sample, value in zip(samples, values):
                writer.writerow([feature, sample, sample_to_habitat[sample], f"{value:.10g}", f"{math.log10(value + 1):.10g}"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="merged_tpm_matrix.txt")
    ap.add_argument("--outdir", default="derepMAG_abundance_stats")
    ap.add_argument("--sample-map", default=None)
    ap.add_argument("--feature-map", default=None, help="Optional TSV with MAG and Feature columns, e.g. K08356 branch labels.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    samples, rows = read_matrix(args.matrix)
    generated_map = os.path.join(args.outdir, "sample_habitat_map.tsv")
    write_sample_map(samples, generated_map)
    sample_map_path = args.sample_map or generated_map
    sample_to_habitat = read_sample_map(sample_map_path, samples)

    total = [0.0] * len(samples)
    mag_values = {}
    for mag, vals in rows:
        mag_values[mag] = vals
        for i, v in enumerate(vals):
            total[i] += v

    def to_groups(values, transform=False):
        groups = defaultdict(list)
        for sample, value in zip(samples, values):
            habitat = sample_to_habitat[sample]
            if habitat in HABITAT_ORDER:
                groups[habitat].append(math.log10(value + 1.0) if transform else value)
        return groups

    total_results = {"all_derepMAG_total_TPM": to_groups(total), "all_derepMAG_log10_TPM_plus_1": to_groups(total, True)}
    write_long(os.path.join(args.outdir, "all_derepMAG_abundance_by_sample.tsv"), samples, sample_to_habitat, {"all_derepMAG_total_TPM": total})
    write_kw(os.path.join(args.outdir, "all_derepMAG_kruskal.tsv"), total_results)
    write_pairwise(os.path.join(args.outdir, "all_derepMAG_pairwise_wilcoxon.tsv"), total_results)

    if args.feature_map:
        feature_to_mags = read_feature_map(args.feature_map)
        feature_values = {}
        unique_feature_mags = sorted({mag for mags in feature_to_mags.values() for mag in mags})
        unique_total = [0.0] * len(samples)
        for mag in unique_feature_mags:
            if mag not in mag_values:
                continue
            for i, v in enumerate(mag_values[mag]):
                unique_total[i] += v
        feature_values["K08356_bearing_unique_MAG_total"] = unique_total
        missing = []
        for feature, mags in feature_to_mags.items():
            vals = [0.0] * len(samples)
            for mag in mags:
                if mag not in mag_values:
                    missing.append((feature, mag))
                    continue
                for i, v in enumerate(mag_values[mag]):
                    vals[i] += v
            feature_values[feature] = vals
        results = {}
        for feature, vals in feature_values.items():
            results[f"{feature}_TPM"] = to_groups(vals)
            results[f"{feature}_log10_TPM_plus_1"] = to_groups(vals, True)
        write_long(os.path.join(args.outdir, "feature_abundance_by_sample.tsv"), samples, sample_to_habitat, feature_values)
        write_kw(os.path.join(args.outdir, "feature_kruskal.tsv"), results)
        write_pairwise(os.path.join(args.outdir, "feature_pairwise_wilcoxon.tsv"), results)
        proportion_values = {}
        for feature, vals in feature_values.items():
            if feature == "K08356_bearing_unique_MAG_total":
                continue
            proportion_values[feature] = [
                (value / denom if denom > 0 else 0.0)
                for value, denom in zip(vals, unique_total)
            ]
        proportion_results = {}
        for feature, vals in proportion_values.items():
            proportion_results[f"{feature}_proportion_of_unique_K08356_MAG_total"] = to_groups(vals)
        write_long(os.path.join(args.outdir, "feature_proportion_by_sample.tsv"), samples, sample_to_habitat, proportion_values)
        write_kw(os.path.join(args.outdir, "feature_proportion_kruskal.tsv"), proportion_results)
        write_pairwise(os.path.join(args.outdir, "feature_proportion_pairwise_wilcoxon.tsv"), proportion_results)
        if missing:
            with open(os.path.join(args.outdir, "feature_map_missing_mags.tsv"), "w", newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(["Feature", "MAG"])
                writer.writerows(missing)


if __name__ == "__main__":
    main()
