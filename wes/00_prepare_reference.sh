#!/usr/bin/env bash
set -euo pipefail

# Reference preparation corresponding to the GRCh38/dbSNP resources used by
# the manuscript WES workflow.
#
# Required:
#   REF_DIR  destination directory.

: "${REF_DIR:?Set REF_DIR to the reference directory}"
mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

if [[ ! -s Homo_sapiens_assembly38.fasta ]]; then
    wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
fi
if [[ ! -s Homo_sapiens_assembly38.fasta.fai ]]; then
    samtools faidx Homo_sapiens_assembly38.fasta
fi
if [[ ! -s Homo_sapiens_assembly38.dict ]]; then
    gatk CreateSequenceDictionary -R Homo_sapiens_assembly38.fasta
fi
if [[ ! -s Homo_sapiens_assembly38.fasta.bwt ]]; then
    bwa index Homo_sapiens_assembly38.fasta
fi

# The executed pipeline ultimately used a chromosome-prefixed dbSNP 138 VCF
# named dbsnp_138.hg38.chr.vcf.gz. The original project included a separate
# chromosome-name conversion step; users should provide an equivalent
# GRCh38/dbSNP138 resource that is compatible with the reference contig names.

echo "Reference FASTA prepared. Provide/index:"
echo "  ${REF_DIR}/dbsnp_138.hg38.chr.vcf.gz"
