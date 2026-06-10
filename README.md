```
--------------------------------------------------------
    888      e      Y88b    /        e88~~\  ~~~888~~~
    888     d8b      Y88b  /        d888        888   
    888    /Y88b      Y88b/         8888 __     888   
    888   /  Y88b     /Y88b         8888   |    888   
|   88P  /____Y88b   /  Y88b        Y888   |    888   
 \__8"  /      Y88b /    Y88b        "88__/     888   

                   5baseTAPS 1.0.0
--------------------------------------------------------
```

# 5baseTAPS: Nextflow Pipeline

A Nextflow pipeline for whole-genome cytosine methylation and SNP variant calling from
**Illumina 5-base** and **TAPS Watchmaker** duplex UMI sequencing data.

Built at [The Jackson Laboratory](https://www.jax.org/) Genome Technologies core, extending
[nf-core/fastquorum](https://nf-co.re/fastquorum) with rastair methylation calling and GATK
HaplotypeCaller variant calling.

---

## Overview

Both library types (Illumina 5-base and TAPS Watchmaker) produce reads with the same
base-level methylation signature: 5-methylcytosine (5mC) appears as thymine (C→T), while
unmethylated cytosines are read as C. The pipeline exploits this shared read representation
to run a single analysis path regardless of library chemistry.

### Workflow

The pipeline is structured as three collaborating Nextflow workflows:

```
FASTQUORUM          UMI extraction → alignment → consensus calling (fgbio)
RASTAIR_METHYLSEQ   CpG methylation calling + M-bias QC + methylation controls (rastair)
GATK_VARIANTCALL    Germline SNP/INDEL calling, CpG-masked, scatter-gather (GATK HC)
```

**Step-by-step:**

1. Raw read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Chunked alignment of raw FASTQs for large samples ([seqkit](https://bioinf.shenwei.me/seqkit/) + [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2))
3. UMI extraction and grouping ([fgbio FastqToBam](http://fulcrumgenomics.github.io/fgbio/), [GroupReadsByUmi](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html))
4. Duplex consensus calling and filtering ([fgbio CallDuplexConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html), [FilterConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/FilterConsensusReads.html))
5. Duplex QC metrics ([fgbio CollectDuplexSeqMetrics](http://fulcrumgenomics.github.io/fgbio/tools/latest/CollectDuplexSeqMetrics.html))
6. Final consensus re-alignment ([bwa-mem2](https://github.com/bwa-mem2/bwa-mem2))
7. CpG methylation calling + per-CpG BED/VCF + per-read BED ([rastair](https://github.com/sbludwig/rastair))
8. M-bias QC report ([rastair](https://github.com/sbludwig/rastair))
9. Lambda (negative) and pUC19 (positive) methylation control QC
10. CpG site mask generation for GATK (derived from the rastair call BED)
11. DRAGstr model calibration + scatter-gather SNP/INDEL calling ([GATK HaplotypeCaller](https://gatk.broadinstitute.org/))
12. MultiQC report combining all QC metrics ([MultiQC](https://multiqc.info/))

---

## Requirements

- **Nextflow** ≥ 24.04.2
- **Singularity** (all processes run in containers; Docker or Apptainer also supported)
- **SLURM** (configured for JAX clusters; other executors require config changes)

---

## Input Samplesheet

Prepare a CSV samplesheet with one row per FASTQ pair (multi-lane samples use multiple rows with the same `sample` name):

```csv
sample,fastq_1,fastq_2,read_structure
SAMPLE1,/path/to/SAMPLE1_R1.fastq.gz,/path/to/SAMPLE1_R2.fastq.gz,7M1S+T 7M1S+T
SAMPLE2,/path/to/SAMPLE2_L1_R1.fastq.gz,/path/to/SAMPLE2_L1_R2.fastq.gz,7M1S+T 7M1S+T
SAMPLE2,/path/to/SAMPLE2_L2_R1.fastq.gz,/path/to/SAMPLE2_L2_R2.fastq.gz,7M1S+T 7M1S+T
```

> **Read structure for Illumina 5-base:** `7M1S+T 7M1S+T`
> (7 bp UMI, 1 bp spacer, then template bases — on both R1 and R2)
>
> See the [fgbio read structure docs](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) for other library configurations.

Multi-lane rows with the same `sample` identifier are merged automatically before UMI grouping.

---

## Running the Pipeline

### Minimal run command

```bash
nextflow run TheJacksonLaboratory/5baseTAPS \
    -profile sumner2_singularity \
    --input samplesheet.csv \
    --genome CHM13 \
    --outdir results/
```

### Recommended: run via SLURM head script

Submit a SLURM wrapper script that calls `nextflow run` so the Nextflow process itself is
managed by the scheduler. Nextflow then submits each pipeline task as a separate SLURM job.

```bash
#!/bin/bash
#SBATCH --job-name=nf-5baseTAPS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x-%j.log
#SBATCH --error=logs/%x-%j.log

module load singularity nextflow

nextflow run TheJacksonLaboratory/5baseTAPS \
    -profile sumner2_singularity \
    --input samplesheet.csv \
    --genome CHM13 \
    --outdir results/ \
    -resume
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | — | Path to samplesheet CSV (required) |
| `--genome` | `CHM13` | Reference genome key (`CHM13` or `GRCh38`; pre-indexed on Elion) |
| `--fasta` | — | Path to reference FASTA (overrides `--genome`) |
| `--fasta_fai` | — | Path to FASTA index (`.fai`) |
| `--dict` | — | Path to sequence dictionary (`.dict`) |
| `--bwamem2` | — | Path to bwa-mem2 index directory |
| `--outdir` | — | Output directory (required) |
| `--mode` | `rd` | Pipeline mode: `rd` (research & development) or `ht` (high throughput) |
| `--run_gatk` | `true` | Run GATK HaplotypeCaller variant calling |
| `--align_raw_bam_chunks` | `1` | Number of parallel alignment chunks per sample (set to 5–8 for large samples ≥100M reads for ~1.83× speedup) |
| `--duplex_seq` | `true` | Enable duplex consensus mode (set `false` for single-strand UMI libraries) |
| `--filter_min_reads` | `1 0 0` | Minimum reads to retain a duplex consensus (format: `AB SS DS`) |
| `--filter_max_base_error_rate` | `0.1` | Maximum per-base error rate for consensus filtering |
| `--rastair_rscript_dir` | `null` | Override rastair R scripts directory (uses bundled `assets/rastair_scripts/` by default) |

---

## Output Structure

```
<outdir>/
├── <sample>/
│   ├── bam/               — duplex consensus BAM + index
│   ├── methylation/       — per-CpG BED/VCF, per-read BED, M-bias HTML, methylKit, summaries
│   ├── variants/          — GATK VCFs, CpG site mask
│   └── qc/                — per-sample QC (alignment, coverage, duplex, FastQC, variant)
├── qc/                    — batch TAPS QC report (all samples combined)
├── report/                — MultiQC HTML report
└── pipeline_info/         — Nextflow execution reports and parameter logs
```

See [docs/output.md](docs/output.md) for a complete description of all output files.

---

## Benchmarking

Validated against Illumina DRAGEN 5-base somatic pipeline (v4.4.6) on a matched whole-genome
sample (GT26-01980, 583M reads, GRCh38):

**SNP calling** (PASS SNPs, chr1–22+X+Y, bcftools norm, vcfeval):

| Caller | Precision | Sensitivity | F1 |
|--------|-----------|-------------|-----|
| Rastair (genome-wide) | 91.8% | 97.2% | **94.4%** |
| GATK HC + CpG mask (non-CpG) | 93.3% | 95.0% | **94.1%** |

**5mC methylation concordance** (54.6 M shared CpG sites, reference CpGs, no variant sites):

| Metric | Value |
|--------|-------|
| Pearson R² | **98.1%** |
| MAE | 1.30 pp |
| Mean bias (DRAGEN − rastair) | +0.16 pp |

See [docs/benchmarking_vcfs.md](docs/benchmarking_vcfs.md) for full benchmarking methods and results.

---

## Documentation

| Document | Contents |
|----------|---------|
| [docs/usage.md](docs/usage.md) | Samplesheet format, parameters, running the pipeline |
| [docs/output.md](docs/output.md) | Complete output file reference |
| [docs/interpret_multiqc.md](docs/interpret_multiqc.md) | How to read the MultiQC report |
| [docs/benchmarking_vcfs.md](docs/benchmarking_vcfs.md) | SNP and methylation concordance vs. Illumina DRAGEN |
| [CITATIONS.md](CITATIONS.md) | Tool citations |

---

## Credits

The JAX-GT pipeline was developed by the JAX Genome Technologies bioinformatics team, built on
top of [nf-core/fastquorum](https://nf-co.re/fastquorum) (Nils Homer & Zach Norgaard,
Fulcrum Genomics) and [rastair](https://github.com/sbludwig/rastair) (Etzioni et al., bioRxiv 2026).
The TAPS methylation conversion subworkflow (`subworkflows/nf-core/bam_taps_conversion`) is
adapted from [nf-core/methylseq](https://nf-co.re/methylseq) (Phil Ewels et al.;
doi:[10.5281/zenodo.1343417](https://doi.org/10.5281/zenodo.1343417)).
