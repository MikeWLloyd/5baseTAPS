# 5baseTAPS MultiQC Report: Interpretation Guide

The MultiQC report (`report/multiqc_report.html`) is the first place to review after a pipeline run. It aggregates QC metrics from every step — raw reads through variant calling — into a single HTML file.

This guide describes each section of the report, the metrics shown, expected ranges, and what to investigate when values fall outside those ranges.

---

## Report Layout

The report opens with a **General Statistics** table at the top, followed by custom TAPS sections pinned below it, then native tool sections (FastQC, mosdepth, bcftools). Sections are ordered roughly by workflow step.

| Approx. position | Section |
|-----------------|---------|
| Top | General Statistics |
| ↓ | Mapping Summary |
| ↓ | Duplex UMI Metrics |
| ↓ | Duplex Family Size Distribution |
| ↓ | CpG Methylation Summary |
| ↓ | Read Partitioning |
| ↓ | CpG Methylation: Sample + Controls |
| ↓ | WGS Coverage |
| ↓ | Variant Calling Summary |
| ↓ | FastQC (per-section plots) |
| ↓ | mosdepth (coverage distribution) |
| Bottom | bcftools stats |

> The M-bias plot and duplex QC PDF are supplementary files not included in the MultiQC report. See the [output reference](./output.md) for their descriptions.

---

## General Statistics

The top table gives one row per sample with the most important cross-tool metrics. Columns are drawn from FastQC, samtools flagstat, mosdepth, and bcftools stats.

| Column | Source | Description |
|--------|--------|-------------|
| Seqs | FastQC | Total raw reads in the FASTQ |
| Avg Len | FastQC | Average read length |
| % GC | FastQC | GC content of raw reads |
| % Aligned | samtools flagstat (pre-dedup) | % of primary reads aligned |
| Mean depth | mosdepth | Mean WGS coverage |
| Median depth | mosdepth | Median WGS coverage |
| ≥ 10x | mosdepth | % of bases covered at ≥ 10× |
| ≥ 30x | mosdepth | % of bases covered at ≥ 30× |
| Max coverage | mosdepth | Peak depth (useful for detecting high-coverage outlier regions) |
| Variants | bcftools | Total PASS variants (SNPs + INDELs) |
| SNPs | bcftools | PASS SNPs |
| INDELs | bcftools | PASS INDELs |
| Ti/Tv | bcftools | Transition/Transversion ratio |

---

## Mapping Summary

> **Source:** pre-dedup samtools flagstat — primary reads before UMI deduplication

These counts are from the pre-deduplication flagstat and reflect raw reads mapped to the sample genome. For deduplication fold-reduction, see the Duplex UMI Metrics section.

| Metric | Expected | Investigate if |
|--------|----------|---------------|
| Mapped % | > 95% | < 90%: possible contamination, wrong reference, or library failure |
| Properly paired % | > 90% | Low paired % can indicate library degradation or mis-specified read groups |
| Unmapped % | < 5% | > 10%: investigate with `samtools view -f 4` to check unaligned read sequences |

---

## Duplex UMI Metrics

> **Source:** fgbio CollectDuplexSeqMetrics + pre/post-dedup flagstat

| Metric | Description | Expected range |
|--------|-------------|---------------|
| Raw reads | Total reads before UMI deduplication | — |
| Consensus reads | Reads after UMI deduplication | — |
| Fold reduction | Raw reads ÷ consensus reads | 1.1–3× for WGS at 30–50× depth |
| Duplex yield % | % of UMI families with reads on both strands (DS families) | Low at WGS depth; higher at ultra-deep; see note below |
| Ideal duplex yield % | Expected duplex yield under ideal equal-strand-sampling model | Benchmark for observed yield |
| DS families | Dual-strand families: source molecules with both strands captured (ds_families) | — |
| True duplexes | Molecules with both strands sequenced (ds_duplexes) | — |
| SS families | Single-strand UMI family observations (ss_families) | SS >> DS is normal for WGS; inverts at very deep sequencing |

> **Duplex yield % context:** For whole-genome 5-base sequencing at 30–50× depth, duplex yield is typically < 5%. It is a function of library depth per molecule — duplex yield increases significantly only at ultra-deep targeted or amplicon sequencing. Low duplex yield at WGS depth is expected, not a problem.
>
> **Fold reduction** reflects how many raw reads were collapsed into one consensus molecule. A fold-reduction of 1.2× means 20% of raw reads were duplicates. Very high fold-reduction (> 5×) indicates a low-complexity library.

---

## Duplex Family Size Distribution

> **Source:** fgbio CollectDuplexSeqMetrics

Bar chart showing the fraction of duplex (DS) families at each UMI family size (1, 2, 3, 4+).

| Pattern | Interpretation |
|---------|---------------|
| Dominated by Size 1 (> 90%) | Normal for WGS at 30–50×; most molecules are captured once per strand |
| Growing Size 2–3 | Expected at higher depth; more reads per UMI = better per-molecule error correction |
| High Size 4+ fraction | Very high duplication; low library complexity or very deep sequencing |
| Very flat distribution | Unusual; may indicate UMI grouping issues |

---

## WGS Coverage

> **Source:** mosdepth on `<sample>.cons.filtered.bam`

| Metric | Expected | Investigate if |
|--------|----------|---------------|
| Mean depth | 30–50× for WGS | < 20×: low input or poor library; > 100×: verify not contaminated with spike-in |
| ≥ 1x % | > 95% | < 90%: significant uncovered regions — check adapter content, reference mismatch |
| ≥ 5x % | > 90% for 30× runs | — |
| ≥ 10x % | > 80% for 30× runs | — |
| ≥ 20x % | > 60% for 30× runs | — |
| chrX depth | ~0.5× of autosomes (XY sample) or ~1× (XX sample) | Use to infer sample sex; unexpected values may indicate sample swap |
| chrY depth | ~0.5× of autosomes (XY) or near 0 (XX) | Near-zero chrY with normal chrX = XX; ~equal chrX/chrY = XY |

> **Sex inference:** Compare chrX and chrY depth to mean autosomal depth. For an XX sample, chrX ≈ autosomal depth and chrY ≈ 0. For an XY sample, chrX ≈ 0.5× autosomal and chrY ≈ 0.5× autosomal.

---

## CpG Methylation Summary

> **Source:** rastair `<sample>.methylation_summary.tsv` (ALL row) + spike-in summaries

| Metric | Expected | Investigate if |
|--------|----------|---------------|
| mCpG % (sample) | 70–80% for human WGS | < 60%: incomplete conversion or sample-specific biology; > 90%: reagent issue or wrong sample |
| λ mCpG % | < 2% | > 2%: incomplete TAPS conversion — bisulfite conversion efficiency analog |
| pUC19 mCpG % | > 90% | < 80%: TAPS reagent failure; pUC19 DNA should be fully methylated |
| CpGs ≥ 1x | ~55–60M for human WGS at 30–50× | Significantly fewer = low coverage; significantly more than genome CpG count = check reference |
| CpGs ≥ 5x % | > 70% for 30× runs | — |
| CpGs ≥ 10x % | > 50% for 30× runs | — |
| CpGs ≥ 20x % | > 20% for 30× runs | — |

> **Lambda conversion efficiency** is the primary TAPS assay QC. The formula is: `Conversion efficiency (%) = 100 − lambda mCpG %`. Lambda should be 0% methylated; any signal indicates incomplete TAPS conversion. Values < 2% (i.e., > 98% conversion) are acceptable. Values of 2–5% are borderline and should be noted in the analysis report. Values > 5% indicate significant conversion failure and may affect methylation calls genome-wide.
>
> **pUC19 positive control** confirms the TAPS chemistry is detecting methylation. pUC19 DNA is in vitro CpG-methylated and should read as ~100% methylated. Values dropping below 90% suggest the TAPS oxidation or borane reduction step was suboptimal.

---

## CpG Methylation: Sample + Controls

> **Source:** same as CpG Methylation Summary

Grouped bar chart showing `mCpG %` for each sample's genome, lambda spike-in, and pUC19 spike-in side by side. Use this for quick visual comparison across multiple samples in a run — the lambda and pUC19 bars should be consistent across samples.

| Bar | Color | Expected |
|-----|-------|----------|
| Sample genome | Blue | 70–80% (human) |
| Lambda (neg ctrl) | Orange | < 2% |
| pUC19 (pos ctrl) | Green | > 90% |

---

## Read Partitioning

> **Source:** pre-dedup samtools flagstat, combining all three genome alignments

Shows what fraction of total **raw** reads were assigned to each genome.

| Category | Expected | Investigate if |
|----------|----------|---------------|
| Sample genome % | 93–98% | < 90%: high spike-in ratio or sample quality issue |
| Lambda (neg ctrl) % | 0.5–2% | > 5%: spike-in over-added during library prep |
| pUC19 (pos ctrl) % | 0.5–3% | > 5%: spike-in over-added |
| Unmapped % | < 5% | > 10%: contamination or wrong reference genome |

> The exact spike-in percentages depend on the amount added during library preparation and can vary by protocol version. What matters is consistency across samples within a run, and that lambda mCpG % remains below 2%.

---

## Variant Calling Summary

> **Source:** bcftools stats on `<sample>.haplotypecaller.filtered-SNP.vcf.gz` (PASS SNPs only, CpG-masked)

| Metric | Expected | Investigate if |
|--------|----------|---------------|
| SNPs | 2.5–3.5M for human WGS (non-CpG only, CpG-masked) | < 1M: severe under-calling; > 5M: likely false positives |
| INDELs | 0.3–0.5M for human WGS | High INDEL count can indicate alignment artifacts |
| Ti/Tv ratio | 2.0–2.1 for WGS germline | < 1.8: elevated false-positive rate (transversions); > 2.5: CpG enrichment artifact or filtering issue |
| Multiallelic sites | < 1% of total variants | High multiallelic % can indicate alignment artifacts at repetitive regions |

> **Ti/Tv context for TAPS:** Because CpG sites are excluded by the CpG mask before GATK runs, the remaining SNP calls are non-CpG variants. CpG transitions (C→T at CpG) normally inflate Ti/Tv in unmasked calls. With the CpG mask applied, Ti/Tv should be close to the genome-wide background of ~2.0–2.1. Values above 2.5 in CpG-masked calls may indicate the mask is not fully removing all TAPS artifacts.
>
> **Comparison to Illumina DRAGEN:** The 5baseTAPS pipeline has been benchmarked at ~94% F1 vs. DRAGEN 5-base on a >580M-read whole-genome sample. See [benchmarking_vcfs.md](./benchmarking_vcfs.md) for the full concordance table.

---

## FastQC Sections

> **Source:** FastQC on raw FASTQ files (pre-alignment)

FastQC is run on raw reads before any processing. Several modules look unusual in TAPS data by design — do not interpret them the same way as standard WGS.

### Per-Base Sequence Quality
Expected: Q30+ across the full read length. A drop at read ends is normal. If quality drops sharply in the first 10–15 bases, check if adapter trimming is needed.

### Per-Base Sequence Content
**Expected to look abnormal in TAPS data.** 5mC→T conversion causes an elevated T fraction at positions corresponding to methylated CpG contexts. This is not an error — it is the expected TAPS read signature. A flat base-content curve across the entire read would be unexpected for a methylated sample.

### Sequence Duplication Levels
PCR duplication is expected in 5-base libraries. The raw duplication level reported by FastQC is based on sequence identity alone; UMI-aware deduplication (fgbio) handles actual duplicates at the molecule level. The **Fold reduction** in the Duplex UMI Metrics section is the relevant duplication metric for this pipeline, not the FastQC duplication level.

### Adapter Content
Should be low (< 5%). Significant adapter contamination indicates short insert sizes; this can be caused by degraded input DNA.

### Per-Sequence GC Content
Expected to show a slightly bimodal or shifted distribution vs. the theoretical genome GC, due to the T-enrichment from TAPS conversion. A sharp narrow spike (possible contamination) or very broad flat distribution (library complexity issue) would warrant investigation.

---

## Mosdepth Coverage Plot

> **Source:** mosdepth on `<sample>.cons.filtered.bam`

The mosdepth section shows per-chromosome coverage distribution plots. Key things to check:

- **Coverage uniformity:** coverage should be roughly uniform across autosomes. Sharp dips at specific chromosomes may indicate regions of low mappability or structural variants.
- **chrX and chrY relative depth:** see WGS Coverage section above for sex inference.
- **Spike in very high coverage regions:** visible as a tail in the distribution; expected to be small for WGS.

---

## bcftools Stats Sections

> **Source:** bcftools stats on `<sample>.haplotypecaller.filtered-SNP.vcf.gz`

The bcftools stats section includes several sub-plots:

| Sub-plot | What to check |
|----------|--------------|
| Variant Quality Distribution | PASS variants should cluster at high quality; a flat distribution or spike at low quality suggests the hard filter cutoffs need review |
| Indel Length Distribution | Small indels (±1–3 bp) should dominate; long indels are rare in germline short-read WGS |
| Substitution Types | C>T and G>A transitions should dominate (CpG-depleted but still the most common transition type). A spike in C>A or A>T transversions may indicate oxidative damage or sequencing artifacts |
| Per-Chromosome SNP Density | Roughly uniform across chromosomes; reduced density on chrX/chrY is normal |

---

*For pipeline output file descriptions and interpretation of supplementary files (M-bias plot, duplex QC PDF), see the [output reference](./output.md).*
