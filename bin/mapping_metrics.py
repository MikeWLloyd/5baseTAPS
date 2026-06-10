#!/usr/bin/env python3
"""
TAPS mapping QC metrics — MultiQC custom_content TSV format.

Combines flagstat files from the main genome, lambda (negative control),
and pUC19 (positive control) to report mapping rates and spike-in composition.

Usage:
    mapping_metrics.py <sample_id> <main.flagstat> <lambda.flagstat> <puc19.flagstat>

Output (to stdout): MultiQC custom_content TSV with comment-header metadata.
  Two sections are written to separate files:
    {sample_id}_read_partitioning_mqc.tsv  — stacked bargraph (spike-in breakdown)
    {sample_id}_mapping_summary_mqc.tsv    — table (mapping rate, properly paired, etc.)
"""

import sys
import re
import csv
import os


def parse_flagstat(path):
    """Parse a samtools flagstat file and return key counts as a dict."""
    result = {"total": 0, "mapped": 0, "properly_paired": 0}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            m = re.match(r"^(\d+)", line)
            if not m:
                continue
            count = int(m.group(1))
            if "total" in line:
                result["total"] = count
            elif "mapped (" in line or (line.endswith("mapped") and "mapped" in line):
                if "mapped" in line and "mate mapped" not in line:
                    result["mapped"] = count
            elif "properly paired" in line:
                result["properly_paired"] = count
    return result

def parse_flagstat_prededup(path):
    """Parse a samtools flagstat fastq2bam prededup file and return key counts as a dict."""
    result = {"total": 0, "mapped": 0, "properly_paired": 0}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            m = re.match(r"^(\d+)", line)
            if not m:
                continue
            count = int(m.group(1))
            if line.endswith("primary"):
                result["total"] = count
            elif "primary mapped" in line:
                result["mapped"] = count
            elif "properly paired" in line:
                result["properly_paired"] = count
    return result


def pct(numerator, denominator, decimals=2):
    """Return a rounded percentage float, or empty string if denominator is zero."""
    try:
        n, d = float(numerator), float(denominator)
        if d == 0:
            return ""
        return round(100.0 * n / d, decimals)
    except (ValueError, TypeError):
        return ""


def write_tsv(path, comment_header, column_headers, data_row):
    """Write a MultiQC custom_content TSV with YAML comment header."""
    with open(path, "w", newline="") as fh:
        for line in comment_header:
            fh.write(f"# {line}\n")
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(column_headers)
        writer.writerow(data_row)


def main():
    if len(sys.argv) != 5:
        sys.exit(
            f"Usage: {sys.argv[0]} <sample_id> <main.flagstat> <lambda.flagstat> <puc19.flagstat>"
        )

    sample_id   = sys.argv[1]
    main_path   = sys.argv[2]
    lambda_path = sys.argv[3]
    puc19_path  = sys.argv[4]

    main_fs   = parse_flagstat_prededup(main_path)
    lambda_fs = parse_flagstat_prededup(lambda_path)
    puc19_fs  = parse_flagstat_prededup(puc19_path)

    total_main            = main_fs["total"]
    mapped_main           = main_fs["mapped"]
    properly_paired_main  = main_fs["properly_paired"]
    unmapped_main         = total_main - mapped_main
    mapped_lambda         = lambda_fs["mapped"]
    mapped_puc19          = puc19_fs["mapped"]

    # Spike-in composition denominator: all reads assigned to any genome
    composition_total = total_main 

    sample_pct   = pct(mapped_main,   composition_total)
    lambda_pct   = pct(mapped_lambda, composition_total)
    puc19_pct    = pct(mapped_puc19,  composition_total)
    unmapped_pct = pct(unmapped_main, composition_total)

    # ── Section 1: Read Partitioning (stacked bargraph) ───────────────────────
    partitioning_header = [
        "id: 'read_partitioning'",
        "section_name: 'Read Partitioning'",
        "description: >",
        "  Fraction of total reads assigned to each genome. Lambda reads form the",
        "  unmethylated negative control (target: < 2%). pUC19 reads form the",
        "  fully-methylated positive control (typically 1–5% of total).",
        "  Unmapped: reads that did not align to any of the three references.",
        "help_text: >",
        "  **Lambda (negative control):** spike-in of unmethylated lambda phage DNA.",
        "  Expected fraction < 2%. High lambda % indicates spike-in over-representation",
        "  or incomplete removal during library prep.",
        "",
        "  **pUC19 (positive control):** spike-in of fully CpG-methylated plasmid.",
        "  Expected 1–5% depending on the spike-in ratio used.",
        "",
        "  **Sample:** all reads mapping to the target genome (main reference).",
        "",
        "  **Unmapped:** reads that did not align to any of the three references.",
        "plot_type: 'bargraph'",
        "pconfig:",
        "  id: 'read_partitioning_plot'",
        "  title: 'Read Partitioning by Genome'",
        "  ylab: '% of total reads'",
        "  cpswitch: false",
        "  stacking: 'normal'",
        "cats:",
        "  Sample genome:",
        "    color: '#2196F3'",
        "  Lambda (neg ctrl):",
        "    color: '#FF9800'",
        "  pUC19 (pos ctrl):",
        "    color: '#4CAF50'",
        "  Unmapped:",
        "    color: '#9E9E9E'",
    ]
    write_tsv(
        f"{sample_id}_read_partitioning_mqc.tsv",
        partitioning_header,
        ["Sample", "Sample genome", "Lambda (neg ctrl)", "pUC19 (pos ctrl)", "Unmapped"],
        [sample_id, sample_pct, lambda_pct, puc19_pct, unmapped_pct],
    )

    # ── Section 2: Mapping Summary (table) ────────────────────────────────────
    mapped_pct_main       = pct(mapped_main,          total_main)
    properly_paired_pct   = pct(properly_paired_main, total_main)
    unmapped_pct_main     = pct(unmapped_main,         total_main)

    summary_header = [
        "id: 'mapping_summary'",
        "section_name: 'Mapping Summary'",
        "description: 'Per-sample read counts and mapping rates against the main reference genome.'",
        "plot_type: 'table'",
        "pconfig:",
        "  id: 'mapping_summary_table'",
        "  title: 'Mapping Summary'",
        "col_config:",
        "  Total reads:",
        "    title: 'Total reads'",
        "    description: 'Total QC-passed reads (main genome flagstat)'",
        "    format: '{:,d}'",
        "  Mapped reads:",
        "    title: 'Mapped'",
        "    description: 'Reads mapped to the main reference genome'",
        "    format: '{:,d}'",
        "  Mapped %:",
        "    title: 'Mapped %'",
        "    description: 'Percentage of total reads mapped to main genome'",
        "    suffix: '%'",
        "    min: 80",
        "    max: 100",
        "    scale: 'RdYlGn'",
        "  Properly paired %:",
        "    title: 'Paired %'",
        "    description: 'Percentage of reads in properly paired alignments'",
        "    suffix: '%'",
        "    scale: 'Blues'",
        "  Unmapped %:",
        "    title: 'Unmapped %'",
        "    description: 'Percentage of reads that did not align to any reference'",
        "    suffix: '%'",
        "    max: 20",
        "    scale: 'RdYlGn_r'",
    ]
    write_tsv(
        f"{sample_id}_mapping_summary_mqc.tsv",
        summary_header,
        ["Sample", "Total R1+R2", "Mapped R1+R2", "Mapped %", "Properly paired %", "Unmapped %"],
        [sample_id, total_main, mapped_main, mapped_pct_main, properly_paired_pct, unmapped_pct_main],
    )


if __name__ == "__main__":
    main()
