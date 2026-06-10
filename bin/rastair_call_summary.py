#!/usr/bin/env python3
"""
Global CpG methylation summary from rastair call BED output.

Two input modes selected by --all-cpgs:

  Default (covered sites only):
    Standard 'rastair call' BED — only positions with at least one read.
    Reports covered CpG count and methylation % over covered sites.
    Columns: sample  chr  total_cpg_sites  global_meth_pct  cpg_cov_ge{1,5,10,20}x

  --all-cpgs (genome-wide):
    BED generated with 'rastair call --cpgs-only --bed-include-empty --all'.
    Columns: sample  chr  global_meth_pct  cpg_cov_ge{1,5,10,20}x

    NOTE: total_cpg_in_ref / covered_cpg_sites / cpg_breadth_pct are omitted.
    The --bed-include-empty BED only covers chromosomes that had reads, making
    total_cpg_in_ref == covered_cpg_sites for partial-genome samples (100% breadth
    is always reported, which is misleading). Drop until a true reference CpG
    count is available.

Usage:
    rastair_call_summary.py <sample_id> <rastair_call.bed.gz> [--per-chr] [--all-cpgs]
"""
# TODO (future): rastair call BED may contain CHH and CHG contexts in addition to CpG.
# Currently this script does not filter by context — it aggregates all sites on the +
# strand, which in practice are CpG-only for the genomes tested. When multi-context
# support is needed, add a context column filter and report CpG / CHH / CHG separately.

import sys
import gzip
import argparse
from collections import defaultdict

COV_THRESHOLDS = [1, 5, 10, 20]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("sample", help="Sample ID")
    p.add_argument("bed_gz", help="rastair_call.bed.gz")
    p.add_argument("--per-chr", action="store_true",
                   help="Also write per-chromosome rows after the global row")
    p.add_argument("--all-cpgs", action="store_true",
                   help="BED was produced with --cpgs-only --bed-include-empty --all; enables genome-wide "
                        "breadth columns (total_cpg_in_ref, covered_cpg_sites, cpg_breadth_pct)")
    return p.parse_args()


def empty_stats():
    return {"mod": 0, "unmod": 0, "total": 0, "cov_counts": defaultdict(int)}


def update(stats, mod, unmod, coverage):
    stats["total"] += 1
    stats["mod"]   += mod
    stats["unmod"] += unmod
    for t in COV_THRESHOLDS:
        if coverage >= t:
            stats["cov_counts"][t] += 1


def methylation_pct(stats):
    denom = stats["mod"] + stats["unmod"]
    return 100.0 * stats["mod"] / denom if denom > 0 else float("nan")


def breadth_pct(stats):
    covered = stats["cov_counts"][1]
    return 100.0 * covered / stats["total"] if stats["total"] > 0 else float("nan")


def print_header(all_cpgs):
    cov_cols = "\t".join(f"cpg_cov_ge{t}x" for t in COV_THRESHOLDS)
    if all_cpgs:
        print(f"sample\tchr\tglobal_meth_pct\t{cov_cols}")
    else:
        print(f"sample\tchr\ttotal_cpg_sites\tglobal_meth_pct\t{cov_cols}")


def print_row(sample, chrom, stats, all_cpgs):
    meth     = methylation_pct(stats)
    cov_vals = "\t".join(str(stats["cov_counts"][t]) for t in COV_THRESHOLDS)
    if all_cpgs:
        print(f"{sample}\t{chrom}\t{meth:.2f}\t{cov_vals}")
    else:
        print(f"{sample}\t{chrom}\t{stats['total']}\t{meth:.2f}\t{cov_vals}")


def main():
    args = parse_args()

    glob = empty_stats()
    per_chr = defaultdict(empty_stats)

    opener = gzip.open if args.bed_gz.endswith(".gz") else open
    with opener(args.bed_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            # columns: chr start end name beta_est strand unmod mod no_snp snp coverage ...
            strand = fields[5]
            if strand != "+":
                continue
            chrom   = fields[0]
            unmod   = int(fields[6])
            mod     = int(fields[7])
            coverage = int(fields[10])

            update(glob, mod, unmod, coverage)
            if args.per_chr:
                update(per_chr[chrom], mod, unmod, coverage)

    print_header(args.all_cpgs)
    print_row(args.sample, "ALL", glob, args.all_cpgs)
    if args.per_chr:
        for chrom in sorted(per_chr):
            print_row(args.sample, chrom, per_chr[chrom], args.all_cpgs)


if __name__ == "__main__":
    main()
