#!/usr/bin/env python3
"""
TAPS methylation QC metrics — MultiQC custom_content TSV format.

Usage:
    taps_methyl_metrics.py <sample_id> <main_summary.tsv> <lambda_summary.tsv> <puc19_summary.tsv>

Outputs:
    {sample_id}_methyl_controls_mqc.tsv  — bargraph: lambda + pUC19 mCpG%
    {sample_id}_methyl_summary_mqc.tsv   — table: main genome methylation + CpG coverage

Pass/fail thresholds:
    lambda_negCtrl : global_meth_pct < 2%   (negative control, should be unmethylated)
    puc19_posCtrl  : global_meth_pct > 90%  (positive control, should be fully methylated)
"""

import sys
import csv


def read_all_row(tsv_path):
    """Return the column dict for the ALL row from a methylation_summary.tsv."""
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("chr") == "ALL":
                return row
    raise ValueError(f"No ALL row found in {tsv_path}")


def pct(value, total, decimals=2):
    try:
        v, t = int(value), int(total)
        if t == 0:
            return ""
        return round(100.0 * v / t, decimals)
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
    if len(sys.argv) < 5:
        sys.exit(
            f"Usage: {sys.argv[0]} <sample_id> <main_summary.tsv> <lambda_summary.tsv> <puc19_summary.tsv> [genome]"
        )

    sample_id  = sys.argv[1]
    main_row   = read_all_row(sys.argv[2])
    lambda_row = read_all_row(sys.argv[3])
    puc19_row  = read_all_row(sys.argv[4])
    genome     = sys.argv[5] if len(sys.argv) > 5 else "Sample"

    lambda_meth = lambda_row.get("global_meth_pct", "")
    puc19_meth  = puc19_row.get("global_meth_pct",  "")
    global_meth = main_row.get("global_meth_pct",   "")

    # Column name for main genome methylation — uses genome name if provided
    genome_col = f"{genome} mCpG %"

    # ── Section 1: CpG methylation grouped bargraph (sample + controls) ─────────
    # 3 side-by-side bars per sample: genome, lambda, pUC19.
    controls_header = [
        "id: 'methyl_controls'",
        "section_name: 'CpG Methylation: Sample + Controls'",
        "description: >",
        "  CpG methylation % for the sample genome and spike-in controls, shown as",
        "  grouped bars (one group per sample, three bars per group).",
        "  Lambda (negative control) target: < 2%. pUC19 (positive control) target: > 90%.",
        "help_text: >",
        f"  **{genome} (sample genome):** CpG methylation % on the reference genome.",
        "  Typical human whole-genome TAPS is 70–80%.",
        "",
        "  **Lambda (negative control):** lambda phage DNA is fully unmethylated.",
        "  Values > 2% indicate incomplete TAPS conversion or reagent contamination.",
        "",
        "  **pUC19 (positive control):** fully CpG-methylated plasmid DNA.",
        "  Values < 90% suggest degraded methylation signal or reagent failure.",
        "plot_type: 'bargraph'",
        "pconfig:",
        "  id: 'methyl_controls_plot'",
        f"  title: 'CpG Methylation: Sample + Controls'",
        "  ylab: 'CpG methylation %'",
        "  ymax: 100",
        "  stacking: 'group'",
        "cats:",
        f"  {genome_col}:",
        "    color: '#2196F3'",
        "  Lambda mCpG %:",
        "    color: '#FF9800'",
        "  pUC19 mCpG %:",
        "    color: '#4CAF50'",
    ]
    write_tsv(
        f"{sample_id}_methyl_controls_mqc.tsv",
        controls_header,
        ["Sample", genome_col, "Lambda mCpG %", "pUC19 mCpG %"],
        [sample_id, global_meth, lambda_meth, puc19_meth],
    )

    # ── Section 2: Main genome methylation + CpG coverage table ──────────────
    global_meth  = main_row.get("global_meth_pct",  "")
    cpg_ge1x     = main_row.get("cpg_cov_ge1x",  "")
    cpg_ge5x     = main_row.get("cpg_cov_ge5x",  "")
    cpg_ge10x    = main_row.get("cpg_cov_ge10x", "")
    cpg_ge20x    = main_row.get("cpg_cov_ge20x", "")

    pct_ge5x  = pct(cpg_ge5x,  cpg_ge1x)
    pct_ge10x = pct(cpg_ge10x, cpg_ge1x)
    pct_ge20x = pct(cpg_ge20x, cpg_ge1x)

    summary_header = [
        "id: 'methyl_summary'",
        "section_name: 'CpG Methylation Summary'",
        "description: 'CpG methylation levels for sample and spike-in controls, plus CpG site coverage breadth.'",
        "plot_type: 'table'",
        "pconfig:",
        "  id: 'methyl_summary_table'",
        "  title: 'CpG Methylation Summary'",
        "col_config:",
        f"  {genome_col}:",
        "    title: 'mCpG %'",
        f"    description: 'CpG methylation % on the {genome} reference genome'",
        "    suffix: '%'",
        "    min: 60",
        "    max: 90",
        "    scale: 'Blues'",
        "  Lambda mCpG %:",
        "    title: 'λ mCpG %'",
        "    description: 'Lambda spike-in methylation % (negative ctrl, target < 2%)'",
        "    suffix: '%'",
        "    max: 2.0",
        "    scale: 'RdYlGn_r'",
        "  pUC19 mCpG %:",
        "    title: 'pUC19 mCpG %'",
        "    description: 'pUC19 spike-in methylation % (positive ctrl, target > 90%)'",
        "    suffix: '%'",
        "    min: 90",
        "    scale: 'RdYlGn'",
        "  CpGs >=1x:",
        "    title: 'CpGs ≥1x'",
        "    description: 'CpG sites with at least 1x coverage (total called)'",
        "    format: '{:,d}'",
        "  PCT CpGs >=5x:",
        "    title: 'CpGs ≥5x %'",
        "    description: '% of covered CpGs (≥1x) with at least 5x coverage'",
        "    suffix: '%'",
        "    scale: 'Blues'",
        "  PCT CpGs >=10x:",
        "    title: 'CpGs ≥10x %'",
        "    description: '% of covered CpGs (≥1x) with at least 10x coverage'",
        "    suffix: '%'",
        "    scale: 'Blues'",
        "  PCT CpGs >=20x:",
        "    title: 'CpGs ≥20x %'",
        "    description: '% of covered CpGs (≥1x) with at least 20x coverage'",
        "    suffix: '%'",
        "    scale: 'Blues'",
    ]
    write_tsv(
        f"{sample_id}_methyl_summary_mqc.tsv",
        summary_header,
        [
            "Sample", genome_col, "Lambda mCpG %", "pUC19 mCpG %",
            "CpGs >=1x", "PCT CpGs >=5x", "PCT CpGs >=10x", "PCT CpGs >=20x",
        ],
        [sample_id, global_meth, lambda_meth, puc19_meth, cpg_ge1x, pct_ge5x, pct_ge10x, pct_ge20x],
    )


if __name__ == "__main__":
    main()
