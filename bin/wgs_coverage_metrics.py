#!/usr/bin/env python3
"""
TAPS WGS coverage metrics — MultiQC custom_content TSV format.

Usage:
    taps_wgs_coverage_metrics.py <sample_id> <mosdepth.summary.txt> \\
        <mosdepth.global.dist.txt> <methylation_summary.tsv>

Outputs:
    {sample_id}_coverage_metrics_mqc.tsv     — table: WGS depth + CpG coverage

mosdepth global.dist.txt (quantized) columns:
    chrom  quantile  fraction_of_bases
Quantile ranges from --quantize 0:1:5:10:20:
    0:1     → 0 coverage
    1:5     → 1–4x
    5:10    → 5–9x
    10:20   → 10–19x
    20:inf  → >=20x
"""

import sys
import csv


def parse_mosdepth_summary(path):
    result = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row.get("chrom", "").strip()
            try:
                result[chrom] = float(row.get("mean", 0))
            except (ValueError, TypeError):
                result[chrom] = 0.0
    return result


def parse_mosdepth_dist(path):
    fractions = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, quantile, fraction = parts[0], parts[1], parts[2]
            if chrom == "total":
                try:
                    fractions[quantile] = float(fraction)
                except (ValueError, TypeError):
                    fractions[quantile] = 0.0
    return fractions


def pct_at_least(dist_fracs, n):
    # Standard mosdepth format: each row is already cumulative (fraction of bases
    # at OR ABOVE that depth). Directly look up the fraction for depth n.
    return dist_fracs.get(str(n), 0.0)


def parse_methyl_summary(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("chr") == "ALL":
                return row
    return {}


def fmt_pct(frac, decimals=2):
    try:
        return round(float(frac) * 100, decimals)
    except (ValueError, TypeError):
        return ""


def cpg_pct(count, total, decimals=2):
    try:
        c, t = int(count), int(total)
        if t == 0:
            return ""
        return round(100.0 * c / t, decimals)
    except (ValueError, TypeError):
        return ""


def write_tsv(path, comment_header, column_headers, data_row):
    with open(path, "w", newline="") as fh:
        for line in comment_header:
            fh.write(f"# {line}\n")
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(column_headers)
        writer.writerow(data_row)


def main():
    if len(sys.argv) != 5:
        sys.exit(
            f"Usage: {sys.argv[0]} <sample_id> <mosdepth.summary.txt> "
            f"<mosdepth.global.dist.txt> <methylation_summary.tsv>"
        )

    sample_id       = sys.argv[1]
    cov_summary     = parse_mosdepth_summary(sys.argv[2])
    dist_fracs      = parse_mosdepth_dist(sys.argv[3])
    methyl_row      = parse_methyl_summary(sys.argv[4])

    mean_cov  = cov_summary.get("total", 0.0)
    mean_cov_str = f"{mean_cov:.4f}" if mean_cov else ""

    pct_1x  = fmt_pct(pct_at_least(dist_fracs, 1))
    pct_5x  = fmt_pct(pct_at_least(dist_fracs, 5))
    pct_10x = fmt_pct(pct_at_least(dist_fracs, 10))
    pct_20x = fmt_pct(pct_at_least(dist_fracs, 20))

    chrx = cov_summary.get("chrX", cov_summary.get("X", ""))
    chry = cov_summary.get("chrY", cov_summary.get("Y", ""))
    chrx_str = f"{float(chrx):.4f}" if chrx != "" else ""
    chry_str = f"{float(chry):.4f}" if chry != "" else ""

    global_meth = methyl_row.get("global_meth_pct", "")
    cpg_ge1x    = methyl_row.get("cpg_cov_ge1x",    "")
    cpg_ge5x    = methyl_row.get("cpg_cov_ge5x",    "")
    cpg_ge10x   = methyl_row.get("cpg_cov_ge10x",   "")
    cpg_ge20x   = methyl_row.get("cpg_cov_ge20x",   "")

    coverage_header = [
        "id: 'coverage_metrics'",
        "section_name: 'WGS Coverage'",
        "description: 'Genome-wide coverage depth from mosdepth on the consensus-aligned BAM.'",
        "help_text: >",
        "  **Mean depth** is computed over all bases in the consensus-aligned BAM.",
        "  Typical target for TAPS whole-genome is 30–50x.",
        "",
        "  **PCT >=Nx** values are read from the mosdepth cumulative global distribution.",
        "",
        "  **chrX/chrY depth** can be used to infer sample sex.",
        "plot_type: 'table'",
        "pconfig:",
        "  id: 'coverage_metrics_table'",
        "  title: 'WGS Coverage'",
        "col_config:",
        "  Mean depth:",
        "    title: 'Depth'",
        "    description: 'Mean alignment coverage across all bases (mosdepth)'",
        "    suffix: 'x'",
        "    min: 0",
        "    scale: 'Blues'",
        "  PCT >=1x:",
        "    title: '>=1x %'",
        "    description: '% of genome bases with at least 1x coverage'",
        "    suffix: '%'",
        "    min: 90",
        "    max: 100",
        "    scale: 'RdYlGn'",
        "  PCT >=5x:",
        "    title: '>=5x %'",
        "    description: '% of genome bases with at least 5x coverage'",
        "    suffix: '%'",
        "    scale: 'Blues'",
        "  PCT >=10x:",
        "    title: '>=10x %'",
        "    description: '% of genome bases with at least 10x coverage'",
        "    suffix: '%'",
        "    scale: 'Blues'",
        "  PCT >=20x:",
        "    title: '>=20x %'",
        "    description: '% of genome bases with at least 20x coverage'",
        "    suffix: '%'",
        "    scale: 'Blues'",
        "  chrX depth:",
        "    title: 'chrX'",
        "    description: 'Mean chrX coverage — use for sex inference'",
        "    suffix: 'x'",
        "  chrY depth:",
        "    title: 'chrY'",
        "    description: 'Mean chrY coverage — use for sex inference'",
        "    suffix: 'x'",
    ]
    write_tsv(
        f"{sample_id}_coverage_metrics_mqc.tsv",
        coverage_header,
        ["Sample", "Mean depth", "PCT >=1x", "PCT >=5x", "PCT >=10x", "PCT >=20x", "chrX depth", "chrY depth"],
        [sample_id, mean_cov_str, pct_1x, pct_5x, pct_10x, pct_20x, chrx_str, chry_str],
    )


if __name__ == "__main__":
    main()
