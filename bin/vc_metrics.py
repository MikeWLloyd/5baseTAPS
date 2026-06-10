#!/usr/bin/env python3
"""
TAPS variant calling metrics — MultiQC custom_content TSV format.

Parses bcftools stats output (from PASS-filtered SNP VCF) and extracts
key variant summary statistics.

Usage:
    taps_vc_metrics.py <sample_id> <bcftools_stats.txt>

Output:
    {sample_id}_variant_calling_mqc.tsv  — table: SNP/INDEL counts + Ti/Tv
"""

import sys
import csv


def parse_bcftools_stats(path):
    result = {
        "snps": 0, "mnps": 0, "indels": 0, "others": 0,
        "multiallelic": 0, "ts": 0, "tv": 0, "tstv_ratio": "",
    }
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("SN"):
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                label = parts[2].strip().rstrip(":")
                try:
                    value = int(parts[3].strip())
                except (ValueError, IndexError):
                    continue
                if   "number of SNPs"               in label: result["snps"]         = value
                elif "number of MNPs"               in label: result["mnps"]         = value
                elif "number of indels"             in label: result["indels"]       = value
                elif "number of others"             in label: result["others"]       = value
                elif "number of multiallelic sites" in label: result["multiallelic"] = value
            elif line.startswith("TSTV"):
                parts = line.split("\t")
                if len(parts) < 6 or parts[2].strip() == "ts":
                    continue
                try:
                    result["ts"]         = int(parts[2].strip())
                    result["tv"]         = int(parts[3].strip())
                    result["tstv_ratio"] = parts[4].strip()
                except (ValueError, IndexError):
                    pass
    return result


def write_tsv(path, comment_header, column_headers, data_row):
    with open(path, "w", newline="") as fh:
        for line in comment_header:
            fh.write(f"# {line}\n")
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(column_headers)
        writer.writerow(data_row)


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <sample_id> <bcftools_stats.txt>")

    sample_id = sys.argv[1]
    stats     = parse_bcftools_stats(sys.argv[2])

    snps         = stats["snps"]
    indels       = stats["indels"]
    multiallelic = stats["multiallelic"]
    tstv_ratio   = stats["tstv_ratio"]
    total_vars   = snps + indels

    try:
        snp_pct   = round(100.0 * snps   / total_vars, 2) if total_vars else ""
        indel_pct = round(100.0 * indels / total_vars, 2) if total_vars else ""
    except ZeroDivisionError:
        snp_pct = indel_pct = ""

    variants_header = [
        "id: 'variant_calling'",
        "section_name: 'Variant Calling Summary'",
        "description: >",
        "  GATK HaplotypeCaller germline variant call summary (PASS SNPs and INDELs).",
        "  Computed from bcftools stats on the hard-filtered, PASS-only VCF.",
        "help_text: >",
        "  **SNPs** and **INDELs** are counted after hard-filtering and restricting",
        "  to PASS variants. SNPs are further filtered with the TAPS CpG mask to",
        "  exclude positions where TAPS 5mC→T conversion could be mistaken for a",
        "  C>T single nucleotide variant.",
        "",
        "  **Ti/Tv ratio:** transition/transversion ratio. Expected ~2.0–2.1 for",
        "  whole-genome germline SNPs in a typical human sample. Lower values may",
        "  indicate false-positive calls; higher values can indicate CpG enrichment.",
        "",
        "  **Multiallelic sites:** positions with >2 alleles; typically a small",
        "  fraction of total variants in germline WGS data.",
        "plot_type: 'table'",
        "pconfig:",
        "  id: 'variant_calling_table'",
        "  title: 'Variant Calling Summary'",
        "col_config:",
        "  Total variants:",
        "    title: 'Total'",
        "    description: 'Total PASS variants (SNPs + INDELs)'",
        "    format: '{:,d}'",
        "  SNPs:",
        "    title: 'SNPs'",
        "    description: 'PASS single nucleotide polymorphisms'",
        "    format: '{:,d}'",
        "  SNP %:",
        "    title: 'SNP %'",
        "    description: 'SNPs as % of total variants'",
        "    suffix: '%'",
        "    scale: 'Blues'",
        "  INDELs:",
        "    title: 'INDELs'",
        "    description: 'PASS insertions and deletions'",
        "    format: '{:,d}'",
        "  Multiallelic:",
        "    title: 'Multiallelic'",
        "    description: 'Sites with more than 2 alleles'",
        "    format: '{:,d}'",
        "  Ti/Tv:",
        "    title: 'Ti/Tv'",
        "    description: 'Transition/Transversion ratio (expected ~2.0–2.1 for WGS)'",
        "    min: 1.8",
        "    max: 2.4",
        "    scale: 'RdYlGn'",
    ]
    write_tsv(
        f"{sample_id}_variant_calling_mqc.tsv",
        variants_header,
        ["Sample", "Total variants", "SNPs", "SNP %", "INDELs", "Multiallelic", "Ti/Tv"],
        [sample_id, total_vars, snps, snp_pct, indels, multiallelic, tstv_ratio],
    )


if __name__ == "__main__":
    main()
