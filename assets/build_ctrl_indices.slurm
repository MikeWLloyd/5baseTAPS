#!/usr/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH --job-name=build_ctrl_bwamem2
#SBATCH --output=build_ctrl_bwamem2.slurm.o
#SBATCH --error=build_ctrl_bwamem2.slurm.e
#SBATCH -p gt_compute

# ── User-configurable paths ────────────────────────────────────────────────
# SIF: bwa-mem2 + samtools singularity image
SIF=/gt/research_development/singularity/cache/bwa-mem2_fgbio_samtools-d6fd27126a192efa.sif

# IGENOMES_BASE: root directory where reference genomes are stored
IGENOMES_BASE=/gt/research_development/igenomes

# Source fasta files to copy from
SOURCE_LAMBDA=/gt/data/seqdma/pipelines/QC_PIPELINE/INDEXES/HUMAN/Bismark_Index_WGBS/Negative_Controls/LamdaDNA/lambda.fasta
SOURCE_PUC19=/gt/data/seqdma/pipelines/QC_PIPELINE/INDEXES/HUMAN/Bismark_Index_WGBS/Positive_Controls/pUC19/pUC19.fasta

# Singularity bind mount (adjust if your data lives on a different mount point)
BIND=/gt:/gt
# ──────────────────────────────────────────────────────────────────────────

module use --append /gt/research_development/qifa/elion/modulefiles
module load singularity

# Derived paths — no need to edit below this line
LAMBDA_WGF="${IGENOMES_BASE}/Controls/Lambda/Sequence/WholeGenomeFasta"
LAMBDA_IDX="${IGENOMES_BASE}/Controls/Lambda/Sequence/BWAmem2Index"
LAMBDA_FA="${LAMBDA_WGF}/genome.fa"

PUC19_WGF="${IGENOMES_BASE}/Controls/pUC19/Sequence/WholeGenomeFasta"
PUC19_IDX="${IGENOMES_BASE}/Controls/pUC19/Sequence/BWAmem2Index"
PUC19_FA="${PUC19_WGF}/genome.fa"

mkdir -p "${LAMBDA_WGF}" "${LAMBDA_IDX}"
mkdir -p "${PUC19_WGF}"  "${PUC19_IDX}"

cp "${SOURCE_LAMBDA}" "${LAMBDA_FA}"
cp "${SOURCE_PUC19}"  "${PUC19_FA}"

# ── Lambda ─────────────────────────────────────────────────────────────────
echo "=== [1/6] bwa-mem2 index for lambda → ${LAMBDA_IDX}/genome.* ==="
singularity exec -B "${BIND}" "${SIF}" \
    bwa-mem2 index -p "${LAMBDA_IDX}/genome" "${LAMBDA_FA}"

echo "=== [2/6] samtools faidx for lambda → ${LAMBDA_FA}.fai ==="
singularity exec -B "${BIND}" "${SIF}" \
    samtools faidx "${LAMBDA_FA}"

echo "=== [3/6] samtools dict for lambda → ${LAMBDA_WGF}/genome.dict ==="
singularity exec -B "${BIND}" "${SIF}" \
    samtools dict "${LAMBDA_FA}" -o "${LAMBDA_WGF}/genome.dict"

# ── pUC19 ──────────────────────────────────────────────────────────────────
echo "=== [4/6] bwa-mem2 index for pUC19 → ${PUC19_IDX}/genome.* ==="
singularity exec -B "${BIND}" "${SIF}" \
    bwa-mem2 index -p "${PUC19_IDX}/genome" "${PUC19_FA}"

echo "=== [5/6] samtools faidx for pUC19 → ${PUC19_FA}.fai ==="
singularity exec -B "${BIND}" "${SIF}" \
    samtools faidx "${PUC19_FA}"

echo "=== [6/6] samtools dict for pUC19 → ${PUC19_WGF}/genome.dict ==="
singularity exec -B "${BIND}" "${SIF}" \
    samtools dict "${PUC19_FA}" -o "${PUC19_WGF}/genome.dict"

# ── Verify outputs ─────────────────────────────────────────────────────────
echo ""
echo "=== Lambda WholeGenomeFasta ==="
ls -lh "${LAMBDA_WGF}/"
echo "=== Lambda BWAmem2Index ==="
ls -lh "${LAMBDA_IDX}/"
echo "=== pUC19 WholeGenomeFasta ==="
ls -lh "${PUC19_WGF}/"
echo "=== pUC19 BWAmem2Index ==="
ls -lh "${PUC19_IDX}/"
