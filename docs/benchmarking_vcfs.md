# 5baseTAPS: Benchmarking vs. Illumina DRAGEN 5-base

**Sample:** GT26-01980 (internal sample, no GIAB truth set available)
**Sequencing depth:** ~583M reads
**Reference:** GRCh38/hg38
**DRAGEN version:** 5-base somatic v4.4.6, min-reads=1 (family depth 1.12 ≈ no UMI consensus)

All comparisons are caller-to-caller concordance (5baseTAPS vs. DRAGEN). Because no orthogonal truth set exists for this sample, concordance is evaluated between callers rather than against a ground truth.

---

## SNP Variant Calling

### Methodology

- All VCFs normalised: `bcftools norm -m-any --check-ref w -f hg38.fa -r chr1-22,X,Y`
- PASS variants only, SNPs only (`bcftools view -f PASS -v snps`)
- Concordance evaluated with vcfeval v3.12.1 (RTG Tools) using `--squash-ploidy`
  (germline GATK vs. somatic DRAGEN — genotype representations differ)
- Chromosomes 1–22, X, Y only

### Variant Counts (chr1-22+X+Y, PASS SNPs, normalised)

| Caller | PASS SNPs |
|--------|-----------|
| DRAGEN (unmasked) | 3,862,394 |
| DRAGEN (CpG-masked) | 2,716,936 |
| Rastair | 4,093,757 |
| GATK HC (CpG-masked) | 2,953,599 |
| Rastair CpG-context SNPs only | 1,191,604 |
| GATK + Rastair CpG-SNPs (combined) | 4,145,203 |

Rastair is a genome-wide caller — it uses three internal ML models (CpG-context, de-novo CpG, and non-CpG "other") and covers all SNP types, not just CpG sites. The CpG-context proportion (~29–30%) is nearly identical between rastair and DRAGEN, reflecting the fixed CpG density of the human genome. The ~231K extra rastair calls correspond to de-novo CpGs (a SNP creates a new CpG context not present in the reference), which DRAGEN does not model.

### SNP Concordance Results

| Arm | Calls | vs. Baseline | Precision | Sensitivity | F1 |
|-----|-------|-------------|-----------|-------------|-----|
| **Arm 1** — Rastair (4.09M genome-wide) | vs. DRAGEN unmasked | 91.76% | 97.15% | **94.37%** |
| **Arm 2** — GATK HC + CpG mask (2.95M non-CpG) | vs. DRAGEN CpG-masked | 93.27% | 95.01% | **94.13%** |
| **Arm 3** — GATK + rastair CpG hybrid (4.15M) | vs. DRAGEN unmasked | 90.52% | 97.11% | **93.70%** |
| **Arm 4** — GATK HC CpG-masked (2.95M) | vs. Rastair CpG-masked | 95.91% | 94.16% | **95.02%** |

#### Notes by Arm

**Arm 1 — Rastair vs. DRAGEN (unmasked):**
Both callers are 5-base-aware and internally separate 5mC→T methylation from true C→T SNPs; no external CpG filter is needed for a fair comparison. This is the direct replication of the [rastair paper (Etzioni et al., bioRxiv 2026)](https://www.biorxiv.org/content/10.1101/2026.03.19.712983) primary comparison. The paper reports F1 = 98.9% restricted to GIAB high-confidence regions (NA12878, 45×); our lower F1 (94.4%) is expected when evaluating over the full genome including low-mappability regions, repeats, and centromeres.

**Arm 2 — GATK HC vs. DRAGEN (both CpG-masked):**
GATK has no knowledge of 5mC→T conversion and would call methylated CpG positions as false-positive C→T SNPs. The CpG mask (28,304,361 intervals derived from the sample's own rastair BED) is applied to both VCFs, recovering ~17 F1 points over naive unmasked GATK (77.9% → 94.1%).

**Arm 3 — Hybrid (GATK non-CpG + rastair CpG-SNPs) vs. DRAGEN (unmasked):**
The hybrid underperforms rastair alone (93.7% vs. 94.4%) because GATK QUAL scores and rastair scores use incompatible scales — vcfeval's threshold sweep preferentially retains GATK calls and discards rastair CpG calls at any non-zero threshold. The all-PASS comparison is the only meaningful one for this combined approach.

**Arm 4 — GATK HC (CpG-masked) vs. Rastair (CpG-masked), non-CpG positions only:**
At positions where methylation ambiguity is removed (CpG mask applied to both), GATK and rastair agree on 95% of SNPs — higher than either achieves individually vs. DRAGEN. This confirms both callers are highly concordant at non-CpG variant sites. Source: `vcfeval_gatk_vs_rastair_masked/summary.txt`, job 2348679.

> **Methods:** All VCFs: PASS SNPs only, `bcftools norm -m-any`, chr1–22+X+Y.

### DRAGEN GermlineStatus composition and its effect on FN rates

DRAGEN somatic mode annotates each PASS SNP with a `GermlineStatus` INFO field reflecting how the variant was classified relative to a population germline database (gnomAD/1000G):

| GermlineStatus | DRAGEN PASS SNPs (unmasked) | % |
|---|---|---|
| `Germline_DB` | 3,723,039 | 96.5% |
| `Somatic` | 136,377 | 3.5% |
| `Germline_Proxi` | 2,978 | 0.1% |

The vast majority of DRAGEN PASS calls are `Germline_DB` — common germline variants present in the population database. Only 3.5% are flagged `Somatic` (not found in the germline DB), which for a non-cancer sample without a matched normal reflects rare or novel germline variants rather than true somatic mutations.

Stratifying the vcfeval TP-baseline and FN by GermlineStatus reveals that the concordance gap is disproportionately concentrated in the `Somatic` category:

**Arm 1 — Rastair vs. DRAGEN (unmasked)**

| GermlineStatus | TP | FN | Total | FN rate | FN share |
|---|---:|---:|---:|---:|---:|
| Germline_DB | 3,645,448 | 77,591 | 3,723,039 | 2.1% | 70.4% |
| Somatic | 104,445 | 31,932 | 136,377 | **23.4%** | 29.0% |
| Germline_Proxi | 2,288 | 690 | 2,978 | 23.2% | 0.6% |

**Arm 2 — GATK HC vs. DRAGEN (both CpG-masked)**

| GermlineStatus | TP | FN | Total | FN rate | FN share |
|---|---:|---:|---:|---:|---:|
| Germline_DB | 2,562,149 | 56,130 | 2,618,279 | 2.1% | 72.8% |
| Somatic | 75,968 | 20,431 | 96,399 | **21.2%** | 26.5% |
| Germline_Proxi | 1,697 | 561 | 2,258 | 24.8% | 0.7% |

Key observations:

- **`Germline_DB` FN rate ≈ 2.1%** in both arms — rastair and GATK are in near-perfect agreement with DRAGEN on common germline SNPs.
- **`Somatic` FN rate ≈ 21–23%** — roughly 1 in 4–5 DRAGEN `Somatic` calls is not recovered by rastair or GATK. This category is **~10× enriched** among FNs relative to its share of the baseline.
- Despite `Somatic` comprising only 3.5% of DRAGEN PASS SNPs, it accounts for ~27–29% of all FNs.

**Interpretation — two competing explanations:**

1. **DRAGEN false positives in the `Somatic` tier.** Without a matched normal, DRAGEN cannot truly distinguish somatic from rare germline variants. Variants absent from the population DB are flagged `Somatic` but may be low-confidence rare germline calls or artefacts. Rastair and GATK rejecting them would be *correct*, not a sensitivity gap — these FNs are DRAGEN FPs that vcfeval incorrectly attributes to the caller under test.

2. **Genuine sensitivity gap at rare variants.** DRAGEN's somatic model is tuned to detect low-VAF minority alleles and may recover rare variants that germline callers (rastair, GATK) filter out for lack of prior support. This would represent a real but narrow sensitivity advantage of DRAGEN at the low end of the allele-frequency spectrum.

Both effects are likely present. Since GT26-01980 is a non-cancer sample, the truly somatic fraction of `Somatic`-labeled calls is expected to be negligible. The practical conclusion is that the headline F1 (~94%) is well-supported at the level of common germline SNPs (>96% of the baseline), and the residual discordance is concentrated in a low-confidence DRAGEN-specific tier whose ground truth cannot be assessed without an orthogonal dataset.

### Ti/Tv ratios across concordance categories

Transition/transversion (Ti/Tv) ratio is an orthogonal quality indicator: real germline SNPs have Ti/Tv ≈ 2.0–2.1 (WGS, all sites) driven by the biochemical transition bias at CpG dinucleotides. Artefacts and random errors approach Ti/Tv ≈ 0.5 (4 transition types vs 8 transversion types by chance). Low Ti/Tv in a call set is therefore a marker of enriched noise.

| Category | n | Ts | Tv | Ti/Tv |
|---|---:|---:|---:|---:|
| **Arm 1 — Rastair vs. DRAGEN (unmasked)** | | | | |
| TP-baseline (concordant) | 3,752,181 | 2,533,187 | 1,218,994 | **2.08** |
| FP (rastair-only) | 337,160 | 186,012 | 150,874 | **1.23** |
| FN (DRAGEN-only, missed by rastair) | 110,213 | 78,247 | 31,966 | **2.45** |
| **Arm 2 — GATK HC vs. DRAGEN (CpG-masked)** | | | | |
| TP-baseline (concordant) | 2,639,814 | 1,622,147 | 1,017,667 | **1.59** |
| FP (GATK-only) | 312,909 | 162,676 | 150,233 | **1.08** |
| FN (DRAGEN-only, missed by GATK) | 77,122 | 47,629 | 29,493 | **1.61** |

**Notes on the TP Ti/Tv difference between arms (2.08 vs. 1.59):** CpG sites drive a large fraction of genome-wide transitions (C→T at methylated CpG is the most frequent germline mutation). Arm 2 excludes CpG positions by design — removing this source of transitions drops the expected WGS Ti/Tv from ~2.1 to ~1.6. Both values are exactly in the expected range for their scope, confirming the TPs are genuine germline SNPs.

**FP Ti/Tv (rastair 1.23, GATK 1.08):** Both are substantially below the TP ratio and far below the germline expectation. GATK's FP set approaches Ti/Tv ≈ 1, close to random noise, indicating the GATK-unique calls not recovered by DRAGEN are largely artefactual. Rastair's FPs are slightly higher (1.23), consistent with de-novo CpG calls (a biology-based transition category rastair models but DRAGEN does not) mixed with low-confidence artefacts.

**FN Ti/Tv — the most informative contrast:**

- **Arm 1 FN Ti/Tv = 2.45** — *elevated* above the TP ratio (2.08). This is the opposite of what you would expect from genuine missed SNPs (which should look like TPs at ~2.08). Elevated Ti/Tv in the FN set means the DRAGEN calls rastair missed are *enriched for C→T transitions*. In 5-base sequencing, 5-methylcytosine is converted to thymine (5mC→T), generating C→T read observations at methylated CpGs. DRAGEN's somatic SNP-calling model may be calling a fraction of these methylation-converted C→T positions as true SNPs. Rastair, being methylation-aware, correctly classifies them as methylation signal and does not emit them as variant calls. This provides additional evidence that part of the DRAGEN `Somatic` FN set represents DRAGEN false positives from 5-base methylation signal rather than genuine SNPs rastair failed to detect.

- **Arm 2 FN Ti/Tv = 1.61** — nearly identical to the TP ratio (1.59). GATK's missed calls look statistically indistinguishable from true concordant variants. This indicates GATK's FNs are genuine sensitivity gaps — real SNPs present in DRAGEN that GATK failed to call — rather than DRAGEN artefacts.

**Summary of Ti/Tv evidence:**

| Observation | Implication |
|---|---|
| TP Ti/Tv ≈ 2.08 / 1.59 (arm-appropriate) | TPs are real germline SNPs |
| FP Ti/Tv ≈ 1.23 / 1.08 (well below TPs) | Caller-unique calls enriched for artefacts |
| Arm 1 FN Ti/Tv = 2.45 (above TPs) | DRAGEN `Somatic` FNs enriched for 5mC→T signal → likely DRAGEN FPs |
| Arm 2 FN Ti/Tv = 1.61 (≈ TPs) | GATK FNs are real SNPs — genuine GATK sensitivity gaps |

The elevated FN Ti/Tv in Arm 1 corroborates the GermlineStatus analysis: the rastair FN set is not composed of real missed SNPs but is substantially contaminated by DRAGEN calls that rastair correctly rejects as methylation signal.

---

## 5mC CpG Methylation Concordance

### Methodology

- **5baseTAPS:** `rastair_call.bed.gz` — `beta_est` field; filtered to reference CpGs (`cpg=REF`), no variant-affected sites (`snp=0`), min coverage 1×, chr1–22+X+Y
- **DRAGEN:** `CX_report.txt.gz` — CG context only (excludes CHG/CHH); β = col3 (methylated) / (col3 + col4); min coverage 1×
- Streaming merge-join on sorted files
- Variant-affected positions excluded on the rastair side (`snp>0`) to avoid sites where DRAGEN CX_report inflates β at het C→T SNPs

### Results

Statistics computed over **54,645,317 shared CpG positions** (intersection of rastair and DRAGEN callable sites, chr1–22+X+Y, reference CpGs only, SNP-affected sites excluded). The intersection covers 96.2% of the union (54.6M / 56.8M) — a highly representative sample.

Methylation at each CpG site is expressed as a beta value β ∈ [0, 1], where 0 = fully unmethylated and 1 = fully methylated.

| Metric | Formula | Value |
|--------|---------|-------|
| Pearson R² | r² of (rastair β, DRAGEN β) pairs | **98.08%** |
| RMSE | √mean((DRAGEN β − rastair β)²) | 0.047 (4.7 pp) |
| MAE | mean(\|DRAGEN β − rastair β\|) | 0.013 (1.30 pp) |
| Mean bias | mean(DRAGEN β − rastair β) | +0.0016 (+0.16 pp, DRAGEN higher) |

> **Interpreting MAE:** At the average CpG site, rastair and DRAGEN differ by 1.30 percentage points in their methylation estimate. The low mean is driven by 94.7% of sites agreeing within 5 pp — the overall 1.30 pp average is pulled up by a small tail of outlier sites.

| Agreement band | % of 54.6 M CpG sites |
|----------------|----------------------|
| \|Δβ\| ≤ 0.05 (5 pp) | 94.73% |
| \|Δβ\| ≤ 0.10 (10 pp) | 98.74% |
| \|Δβ\| ≤ 0.20 (20 pp) | 99.49% |

> The small positive bias (+0.16 pp, DRAGEN slightly higher) is consistent with residual uncorrected het C→T sites in the DRAGEN CX report, as reported in the [rastair paper (Etzioni et al., bioRxiv 2026)](https://www.biorxiv.org/content/10.1101/2026.03.19.712983).

---

## Summary

| Comparison | Metric | Value |
|------------|--------|-------|
| GATK vs. DRAGEN, unmasked (naive) | F1 | 77.9% |
| GATK vs. DRAGEN, CpG-masked | F1 | **94.1%** |
| Rastair vs. DRAGEN, all sites | F1 | **94.4%** |
| GATK + Rastair CpG-SNPs vs. DRAGEN | F1 | 93.7% |
| GATK vs. Rastair, CpG-masked (non-CpG only) | F1 | **95.0%** |
| Rastair vs. DRAGEN, methylation R² | Pearson R² | **98.1%** (n=54,645,317) |
| Rastair vs. DRAGEN, methylation MAE | per-site \|Δβ\| | 0.013 (1.30 pp) |
| Rastair vs. DRAGEN, methylation bias | mean Δβ (DRAGEN − rastair) | +0.0016 (+0.16 pp) |

**Key findings:**

1. **CpG masking is essential** for fair 5-base SNP comparison: without it GATK appears to have F1 = 77.9% (false negatives from 5mC→T conversion positions GATK correctly excluded); with matching CpG masks applied to both callers, GATK achieves F1 = **94.1%**.
2. **Rastair F1 = 94.4%** over the full genome (all regions including difficult). The rastair paper reports 98.9% restricted to GIAB high-confidence regions; the gap reflects whole-genome evaluation scope.
3. **Combined GATK + rastair CpG-SNPs (F1 93.7%)** is slightly lower than rastair alone, suggesting rastair's integrated model is marginally cleaner than the GATK/rastair hybrid.
4. **GATK and rastair agree at 95% F1** at non-CpG positions (Arm 4, both CpG-masked), confirming high caller concordance where methylation ambiguity is removed.
5. **Methylation R² = 98.1%** at 54.6M shared reference CpGs confirms near-perfect genome-wide concordance. Mean bias of +0.16 pp (DRAGEN slightly higher) matches the paper's observation of residual uncorrected het C→T sites in DRAGEN's CX_report.
6. **The F1 gap is driven by DRAGEN's `Somatic` tier.** 96.5% of DRAGEN PASS SNPs are `Germline_DB`; rastair and GATK achieve a FN rate of only ~2.1% on these. The remaining 3.5% (`Somatic` — rare/novel variants absent from the germline DB) carry a ~21–23% FN rate and account for ~27–29% of all FNs. Without a matched normal, DRAGEN cannot distinguish true somatic from rare germline; these discordant calls may reflect DRAGEN false positives rather than rastair/GATK sensitivity gaps.

---

**See also — independent benchmarking with GIAB truth set:** The rastair team benchmarks rastair vs. DRAGEN 5-base on NA12878 (GIAB truth set) and reports rastair F1 = 0.906 vs. DRAGEN F1 = 0.899 — rastair slightly outperforms DRAGEN when orthogonal ground truth is available. Our whole-genome caller-to-caller concordance of 94.4% is consistent with this result. See [rastair.com](https://www.rastair.com/) for details.
