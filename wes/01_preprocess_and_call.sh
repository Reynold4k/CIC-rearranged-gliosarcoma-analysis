#!/usr/bin/env bash
set -euo pipefail

# Reproducible WES workflow corresponding to the analysis used in the manuscript.
#
# Required environment variables:
#   RAW_DATA   Directory containing:
#                Blood/Blood.R1.fq.gz
#                Blood/Blood.R2.fq.gz
#                Primary_CIC/P_CIC.R1.fq.gz
#                Primary_CIC/P_CIC.R2.fq.gz
#                Recurrent-CIC/R_CIC.R1.fq.gz
#                Recurrent-CIC/R_CIC.R2.fq.gz
#   WORK_DIR   Analysis working directory.
#
# Optional:
#   REF_DIR    Reference directory (default: ${WORK_DIR}/reference)
#   OUT_DIR    Output directory (default: ${WORK_DIR}/results)
#   THREADS    Mutect2 PairHMM threads (default: 32)
#   THREADS_PER_SAMPLE  Threads for per-sample preprocessing (default: 10)
#
# The workflow intentionally reproduces the simplified analysis used in the
# manuscript. It does NOT add gnomAD, a panel of normals, contamination
# estimation, or additional indel known-sites resources that were not used in
# the executed analysis.

: "${RAW_DATA:?Set RAW_DATA to the directory containing WES FASTQs}"
: "${WORK_DIR:?Set WORK_DIR to the WES working directory}"

REF_DIR="${REF_DIR:-${WORK_DIR}/reference}"
OUT_DIR="${OUT_DIR:-${WORK_DIR}/results}"
THREADS="${THREADS:-32}"
THREADS_PER_SAMPLE="${THREADS_PER_SAMPLE:-10}"

REF_GENOME="${REF_DIR}/Homo_sapiens_assembly38.fasta"
DBSNP="${REF_DIR}/dbsnp_138.hg38.chr.vcf.gz"

NORMAL_NAME="Blood"
TUMOR1_NAME="Primary_CIC"
TUMOR2_NAME="Recurrent_CIC"

mkdir -p "${OUT_DIR}"/{qc,trimmed,aligned,dedup,bqsr,variants,annotated,reports}

for f in "${REF_GENOME}" "${DBSNP}"; do
    [[ -s "${f}" ]] || { echo "Missing required reference: ${f}" >&2; exit 1; }
done

process_sample() {
    local sample="$1"
    local r1="$2"
    local r2="$3"
    local threads="$4"

    [[ -s "${r1}" ]] || { echo "Missing FASTQ: ${r1}" >&2; exit 1; }
    [[ -s "${r2}" ]] || { echo "Missing FASTQ: ${r2}" >&2; exit 1; }

    echo "[$(date)] Processing ${sample}"

    fastqc -t "${threads}" -o "${OUT_DIR}/qc" "${r1}" "${r2}"

    local out_r1="${OUT_DIR}/trimmed/${sample}_R1.trimmed.fq.gz"
    local out_r2="${OUT_DIR}/trimmed/${sample}_R2.trimmed.fq.gz"
    fastp \
        -i "${r1}" -I "${r2}" \
        -o "${out_r1}" -O "${out_r2}" \
        --thread "${threads}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --n_base_limit 5 \
        --length_required 50 \
        --correction \
        --html "${OUT_DIR}/qc/${sample}_fastp.html" \
        --json "${OUT_DIR}/qc/${sample}_fastp.json" \
        2>&1 | tee "${OUT_DIR}/qc/${sample}_fastp.log"

    local rg="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}_lib1"
    local sorted_bam="${OUT_DIR}/aligned/${sample}.sorted.bam"
    bwa mem -t "${threads}" -R "${rg}" "${REF_GENOME}" "${out_r1}" "${out_r2}" | \
        samtools view -@ 4 -Sb - | \
        samtools sort -@ 4 -o "${sorted_bam}" -
    samtools index -@ 4 "${sorted_bam}"
    samtools flagstat "${sorted_bam}" > "${OUT_DIR}/aligned/${sample}.flagstat.txt"

    local dedup_bam="${OUT_DIR}/dedup/${sample}.dedup.bam"
    gatk MarkDuplicates \
        -I "${sorted_bam}" \
        -O "${dedup_bam}" \
        -M "${OUT_DIR}/dedup/${sample}.metrics.txt" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT

    local recal_table="${OUT_DIR}/bqsr/${sample}.recal_data.table"
    local bqsr_bam="${OUT_DIR}/bqsr/${sample}.bqsr.bam"
    gatk BaseRecalibrator \
        -R "${REF_GENOME}" \
        -I "${dedup_bam}" \
        --known-sites "${DBSNP}" \
        -O "${recal_table}"

    gatk ApplyBQSR \
        -R "${REF_GENOME}" \
        -I "${dedup_bam}" \
        --bqsr-recal-file "${recal_table}" \
        -O "${bqsr_bam}"
    samtools index "${bqsr_bam}"
}

process_sample "Blood" \
    "${RAW_DATA}/Blood/Blood.R1.fq.gz" \
    "${RAW_DATA}/Blood/Blood.R2.fq.gz" \
    "${THREADS_PER_SAMPLE}"

process_sample "Primary_CIC" \
    "${RAW_DATA}/Primary_CIC/P_CIC.R1.fq.gz" \
    "${RAW_DATA}/Primary_CIC/P_CIC.R2.fq.gz" \
    "${THREADS_PER_SAMPLE}"

process_sample "Recurrent_CIC" \
    "${RAW_DATA}/Recurrent-CIC/R_CIC.R1.fq.gz" \
    "${RAW_DATA}/Recurrent-CIC/R_CIC.R2.fq.gz" \
    "${THREADS_PER_SAMPLE}"

BLOOD_BQSR="${OUT_DIR}/bqsr/Blood.bqsr.bam"
PRIMARY_BQSR="${OUT_DIR}/bqsr/Primary_CIC.bqsr.bam"
RECURRENT_BQSR="${OUT_DIR}/bqsr/Recurrent_CIC.bqsr.bam"

gatk Mutect2 \
    -R "${REF_GENOME}" \
    -I "${PRIMARY_BQSR}" \
    -I "${BLOOD_BQSR}" \
    -tumor Primary_CIC \
    -normal Blood \
    --native-pair-hmm-threads "${THREADS}" \
    -O "${OUT_DIR}/variants/Primary_CIC_vs_Blood.vcf.gz"

gatk Mutect2 \
    -R "${REF_GENOME}" \
    -I "${RECURRENT_BQSR}" \
    -I "${BLOOD_BQSR}" \
    -tumor Recurrent_CIC \
    -normal Blood \
    --native-pair-hmm-threads "${THREADS}" \
    -O "${OUT_DIR}/variants/Recurrent_CIC_vs_Blood.vcf.gz"

for pair in Primary_CIC_vs_Blood Recurrent_CIC_vs_Blood; do
    gatk FilterMutectCalls \
        -R "${REF_GENOME}" \
        -V "${OUT_DIR}/variants/${pair}.vcf.gz" \
        -O "${OUT_DIR}/variants/${pair}.filtered.vcf.gz"

    bcftools view -f PASS \
        "${OUT_DIR}/variants/${pair}.filtered.vcf.gz" \
        -Oz -o "${OUT_DIR}/variants/${pair}.filtered.PASS.vcf.gz"
    bcftools index -t "${OUT_DIR}/variants/${pair}.filtered.PASS.vcf.gz"
done

PRIMARY_COUNT=$(bcftools view -H "${OUT_DIR}/variants/Primary_CIC_vs_Blood.filtered.PASS.vcf.gz" | wc -l)
RECURRENT_COUNT=$(bcftools view -H "${OUT_DIR}/variants/Recurrent_CIC_vs_Blood.filtered.PASS.vcf.gz" | wc -l)

{
    echo "Primary_CIC PASS variants: ${PRIMARY_COUNT}"
    echo "Recurrent_CIC PASS variants: ${RECURRENT_COUNT}"
} > "${OUT_DIR}/reports/variant_summary.txt"

multiqc "${OUT_DIR}" -o "${OUT_DIR}/reports" -n WES_analysis_report

echo "WES calling complete."
