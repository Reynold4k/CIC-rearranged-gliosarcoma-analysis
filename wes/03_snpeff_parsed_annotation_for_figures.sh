#!/usr/bin/env bash
set -euo pipefail

# Create the parsed SnpEff tables consumed by the publication-figure scripts.
#
# Required:
#   WORK_DIR  WES analysis root containing results/variants.
#
# Optional:
#   OUT_DIR   default: ${WORK_DIR}/results
#   SNPEFF_DB default: GRCh38.99

: "${WORK_DIR:?Set WORK_DIR to the WES working directory}"
OUT_DIR="${OUT_DIR:-${WORK_DIR}/results}"
ANNOT_DIR="${OUT_DIR}/annotated"
SNPEFF_DB="${SNPEFF_DB:-GRCh38.99}"
mkdir -p "${ANNOT_DIR}"

annotate_sample() {
    local pair_name="$1"
    local input_vcf="${OUT_DIR}/variants/${pair_name}.filtered.PASS.vcf.gz"
    local out_vcf="${ANNOT_DIR}/${pair_name}.snpeff.vcf.gz"
    local out_html="${ANNOT_DIR}/${pair_name}.snpeff_summary.html"
    local out_table="${ANNOT_DIR}/${pair_name}.annotated_variants.tsv"

    snpEff -Xmx24g -v "${SNPEFF_DB}" \
        -stats "${out_html}" \
        "${input_vcf}" \
        2> "${ANNOT_DIR}/${pair_name}.snpeff.log" \
        | bgzip > "${out_vcf}"
    tabix -p vcf "${out_vcf}"

    printf 'CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT\tHGVS_C\tHGVS_P\n' > "${out_table}"

    zcat "${out_vcf}" | grep -v '^#' | while IFS=$'\t' read -r chrom pos id ref alt rest; do
        ann=$(echo "${rest}" | grep -oP 'ANN=[^;]+' | head -1 | sed 's/ANN=//')
        first_ann=$(echo "${ann}" | cut -d',' -f1)
        gene=$(echo "${first_ann}" | cut -d'|' -f4)
        effect=$(echo "${first_ann}" | cut -d'|' -f2)
        impact=$(echo "${first_ann}" | cut -d'|' -f3)
        hgvs_c=$(echo "${first_ann}" | cut -d'|' -f10)
        hgvs_p=$(echo "${first_ann}" | cut -d'|' -f11)
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${chrom}" "${pos}" "${ref}" "${alt}" \
            "${gene}" "${effect}" "${impact}" "${hgvs_c}" "${hgvs_p}"
    done >> "${out_table}"
}

annotate_sample "Primary_CIC_vs_Blood"
annotate_sample "Recurrent_CIC_vs_Blood"
