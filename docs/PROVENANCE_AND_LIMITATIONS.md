# Provenance and analysis limitations

## Why this repository contains both cleaned code and provenance files

The original analysis was performed interactively on an HPC system and was distributed across shell scripts, a long R analysis file, saved Seurat objects and final figure scripts. The manuscript-facing scripts in `wes/` and `scrna/` reorganize only the relevant code into readable analysis stages. They do not introduce new biological analyses.

Selected source files are retained under `provenance/` to show how the cleaned code relates to the recovered analysis.

## WES

The main executed WES workflow used:

- FastQC and fastp
- BWA-MEM alignment to GRCh38
- SAMtools
- GATK MarkDuplicates
- GATK BaseRecalibrator / ApplyBQSR
- matched tumor-normal GATK Mutect2
- FilterMutectCalls
- PASS extraction with bcftools
- SnpEff GRCh38.99 annotation
- exact primary/recurrent comparison by `CHROM:POS:REF:ALT`
- AD-based VAF calculation in the figure code

The workflow was a simplified Mutect2 analysis. The executed calling script did not supply gnomAD, a panel of normals or formal contamination estimation, and BQSR used dbSNP alone. These omissions are intentionally not "corrected" in the public code because doing so would no longer reproduce the analysis reported in the manuscript.

Historical environment setup scripts do not perfectly match versions observed in execution logs. See `SOFTWARE_VERSIONS.md`.

## scRNA-seq

The manuscript-specific analysis starts in the original long R file at the section marked `## whx`. Earlier sections of that R file concern unrelated analyses and were intentionally excluded from this repository.

### QC and clustering

The recovered code directly records the QC and clustering parameters now used in the cleaned scripts. No additional doublet removal or upper feature/UMI thresholds were added.

### Tumor and myeloid annotation

Tumor and myeloid states were manually assigned from reclustered Seurat clusters using marker-expression profiles.

A late-stage plotting script refined the **display mapping** of the myeloid UMAP relative to the earlier mapping used when the downstream refined object was assembled. In particular, the recovered files contain a mapping difference for myeloid clusters 2 and 11. This repository keeps that distinction explicit:

- `02_tumor_myeloid_reclustering.R` reflects the mapping present in the original downstream analysis object used by the recovered pseudotime/CellChat code.
- `06_final_umap_redraw.R` reproduces the later final UMAP display mapping.

This is documented rather than silently harmonized because changing the earlier labels retrospectively would alter the provenance of downstream analyses.

### Pseudotime

Slingshot inferred more than one lineage in the recovered workflow. The manuscript figure/Methods discuss the first lineage; the cleaned script therefore stores the inferred pseudotime matrix but uses lineage 1 for the reported UMAP and Early/Middle/Late heatmap.

Pseudotime is an exploratory transcriptional ordering and is not evidence of chronological or clonal evolution.

### CellChat

Two CellChat analyses were present in the original workflow:

1. all refined cell types;
2. a recomputed network after excluding T/NK and dendritic cells.

The focused ligand-receptor and pathway plots were derived from the first analysis, whereas the tumor-to-myeloid and reciprocal myeloid-to-tumor strength matrices were derived from the second analysis. The repository preserves this distinction.

CellChat predictions are transcript-expression-based computational inferences and should not be interpreted as direct evidence of functional signaling.

## Data and privacy

No raw FASTQs, BAMs, VCFs, large RDS files, microscopy images, clinical identifiers or patient-level source records are included. Absolute HPC paths and the original email notification address were removed from the manuscript-facing scripts.

### Marker-heatmap annotation colors

The recovered manuscript-specific R section referenced a `module_col_map` object that was not defined within the recovered WHX code block. This affects only the small module-annotation color strip, not the averaged expression matrix, z-score calculation, marker definitions or cell annotations. The cleaned script therefore assigns deterministic display-only colors while preserving the analytical values.

## Public-release sanitization

The provenance copies in this repository are intentionally sanitized. Institution-internal HPC paths, allocation identifiers, usernames, and internal sample identifiers were replaced with generic placeholders. No computational settings or analysis decisions were changed during sanitization.
