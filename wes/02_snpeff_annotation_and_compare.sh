#!/usr/bin/env bash
set -euo pipefail

# SnpEff annotation and exact primary-versus-recurrent comparison.
#
# Required:
#   WORK_DIR  WES analysis root containing results/variants.
#
# Optional:
#   OUT_DIR   default: ${WORK_DIR}/results
#   SNPEFF_DB default: GRCh38.99

: "${WORK_DIR:?Set WORK_DIR to the WES working directory}"
OUT_DIR="${OUT_DIR:-${WORK_DIR}/results}"
VARIANT_DIR="${OUT_DIR}/variants"
ANNOT_DIR="${OUT_DIR}/annotated"
REPORT_DIR="${OUT_DIR}/reports"
SNPEFF_DB="${SNPEFF_DB:-GRCh38.99}"

mkdir -p "${ANNOT_DIR}" "${REPORT_DIR}"

PRIMARY_VCF="${VARIANT_DIR}/Primary_CIC_vs_Blood.filtered.PASS.vcf.gz"
RECURRENT_VCF="${VARIANT_DIR}/Recurrent_CIC_vs_Blood.filtered.PASS.vcf.gz"

for f in "${PRIMARY_VCF}" "${RECURRENT_VCF}"; do
    [[ -s "${f}" ]] || { echo "Missing PASS VCF: ${f}" >&2; exit 1; }
done

snpEff download -v "${SNPEFF_DB}" || true

annotate() {
    local label="$1"
    local vcf="$2"
    snpEff -Xmx32g "${SNPEFF_DB}" \
        -stats "${ANNOT_DIR}/${label}.snpEff_summary.html" \
        -csvStats "${ANNOT_DIR}/${label}.snpEff_summary.csv" \
        "${vcf}" | bgzip > "${ANNOT_DIR}/${label}.annotated.vcf.gz"
    bcftools index -t "${ANNOT_DIR}/${label}.annotated.vcf.gz"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/ANN\n' \
        "${ANNOT_DIR}/${label}.annotated.vcf.gz" \
        > "${ANNOT_DIR}/${label}.variants.tsv"
}

annotate "Primary_CIC" "${PRIMARY_VCF}"
annotate "Recurrent_CIC" "${RECURRENT_VCF}"

bcftools query -f '%CHROM:%POS:%REF:%ALT\n' \
    "${ANNOT_DIR}/Primary_CIC.annotated.vcf.gz" | sort -u \
    > "${ANNOT_DIR}/primary_sites.txt"
bcftools query -f '%CHROM:%POS:%REF:%ALT\n' \
    "${ANNOT_DIR}/Recurrent_CIC.annotated.vcf.gz" | sort -u \
    > "${ANNOT_DIR}/recurrent_sites.txt"

comm -12 "${ANNOT_DIR}/primary_sites.txt" "${ANNOT_DIR}/recurrent_sites.txt" \
    > "${ANNOT_DIR}/shared_variants.txt"
comm -23 "${ANNOT_DIR}/primary_sites.txt" "${ANNOT_DIR}/recurrent_sites.txt" \
    > "${ANNOT_DIR}/primary_specific_variants.txt"
comm -13 "${ANNOT_DIR}/primary_sites.txt" "${ANNOT_DIR}/recurrent_sites.txt" \
    > "${ANNOT_DIR}/recurrent_specific_variants.txt"

{
    echo "Primary PASS exact sites: $(wc -l < "${ANNOT_DIR}/primary_sites.txt")"
    echo "Recurrent PASS exact sites: $(wc -l < "${ANNOT_DIR}/recurrent_sites.txt")"
    echo "Exact shared sites: $(wc -l < "${ANNOT_DIR}/shared_variants.txt")"
    echo "Primary-specific sites: $(wc -l < "${ANNOT_DIR}/primary_specific_variants.txt")"
    echo "Recurrent-specific sites: $(wc -l < "${ANNOT_DIR}/recurrent_specific_variants.txt")"
} | tee "${REPORT_DIR}/variant_comparison_summary.txt"
