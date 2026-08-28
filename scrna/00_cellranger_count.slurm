#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G

# Cell Ranger preprocessing used for the recurrent tumor scRNA-seq dataset.
#
# Set these variables before submission:
#   CELLRANGER_REF=/path/to/refdata-gex-GRCh38-2024-A
#   FASTQ_DIR=/path/to/fastqs
# Optional:
#   SAMPLE_ID=RECURRENT_TUMOR
#   OUT_ID=recurrent_tumor_cellranger
#
# The original analysis used Cell Ranger 8.0.1, 24 local cores,
# 128 GB local memory, and --create-bam=true.

module purge
module load CellRanger/8.0.1

: "${CELLRANGER_REF:?Set CELLRANGER_REF}"
: "${FASTQ_DIR:?Set FASTQ_DIR}"
SAMPLE_ID="${SAMPLE_ID:-RECURRENT_TUMOR}"
OUT_ID="${OUT_ID:-recurrent_tumor_cellranger}"

cellranger count \
  --id="${OUT_ID}" \
  --transcriptome="${CELLRANGER_REF}" \
  --fastqs="${FASTQ_DIR}" \
  --sample="${SAMPLE_ID}" \
  --localcores=24 \
  --localmem=128 \
  --create-bam=true
