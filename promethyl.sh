#!/bin/bash
#SBATCH --job-name=promethyl
#SBATCH --partition=prod_short
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=02:00:00
#SBATCH --mem=8GB

# load modules and environment
module load miniconda3
conda activate promethyl

# Paths inside repo
MODKIT_DIR="output/modkit"
ANNOTATION="/home/tim.green/promethyl/reference/CpGs_with_promoters.bed"
INCLUDE_BED="/home/tim.green/promethyl/reference/CPGIslandsBED3.bed"
OUTPUT="output/probands_modkit_cohort_methylation.tsv"
CONFIG="cohort.yaml"

# External reference genome (NOT in repo)
REF="/hpc/bpipeLibrary/shared/cpipe-wgs/hg38/v0/Homo_sapiens_assembly38.fasta"

# Threads (fixed for Option A)
THREADS=8

echo "==================================="
echo "Promethyl pipeline starting"
echo "Config: $CONFIG"
echo "Repo root: $ROOT_DIR"
echo "Threads: $THREADS"
echo "Start: $(date)
echo "==================================="

START=$(date +%s)

# Ensure output directories exist
mkdir -p "output/modkit"

# Run pipeline
python "/home/tim.green/promethyl/src/main.py" \
    --config "$CONFIG" \
    --modkit-dir "$MODKIT_DIR" \
    --annotation "$ANNOTATION" \
    --ref "$REF" \
    --output "$OUTPUT" \
    --include-bed "$INCLUDE_BED" \
    --threads "$THREADS"

END=$(date +%s)
ELAPSED=$(( END - START ))
HOURS=$(( ELAPSED / 3600 ))
MINUTES=$(( (ELAPSED % 3600) / 60 ))
SECONDS=$(( ELAPSED % 60 ))

echo "==================================="
echo "Pipeline finished successfully"
echo "End     : $(date)"
echo "Runtime : ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo "==================================="