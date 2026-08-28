# CIC-rearranged gliosarcoma analysis code

This repository contains the analysis code used for the longitudinal whole-exome sequencing (WES) and single-cell RNA-sequencing (scRNA-seq) analyses in a manuscript describing a primary and recurrent CIC-rearranged gliosarcoma.

The repository was assembled from the original HPC analysis scripts and the final figure-generation code. The aim is to document **what was actually done**, rather than to replace the original workflow with a different "best-practice" pipeline.

## Scope

The repository covers the computational analyses that are represented in the manuscript:

- WES read QC, trimming, GRCh38 alignment, duplicate marking and base-quality score recalibration
- matched tumor-normal somatic variant calling with GATK Mutect2
- PASS-only somatic variant filtering
- SnpEff functional annotation
- exact longitudinal comparison using `CHROM + POS + REF + ALT`
- AD-based variant allele fraction (VAF) calculations and WES publication figures
- Cell Ranger processing of the recurrent-tumor scRNA-seq data
- Seurat QC, normalization, PCA, graph clustering and UMAP
- major cell-type annotation
- tumor-cell and macrophage/myeloid reclustering
- marker-based annotation and representative marker heatmaps
- relative myeloid-state composition
- exploratory Slingshot pseudotime analysis
- CellChat ligand-receptor and pathway inference
- final UMAP rendering used for manuscript figures

Multiplex immunofluorescence is an experimental imaging assay and is therefore not included in this computational code repository.

## Repository structure

```text
.
├── README.md
├── docs/
│   ├── ANALYSIS_TO_MANUSCRIPT_MAP.md
│   ├── METHODS_EVIDENCE_SUMMARY.md
│   ├── PROVENANCE_AND_LIMITATIONS.md
│   └── SOFTWARE_VERSIONS.md
├── wes/
│   ├── 00_prepare_reference.sh
│   ├── 01_preprocess_and_call.sh
│   ├── 02_snpeff_annotation_and_compare.sh
│   ├── 03_snpeff_parsed_annotation_for_figures.sh
│   ├── generate_wes_publication_figures.py
│   └── generate_exact_shared_vaf_dumbbell.py
├── scrna/
│   ├── 00_cellranger_count.slurm
│   ├── 01_global_seurat_qc_clustering.R
│   ├── 02_tumor_myeloid_reclustering.R
│   ├── 03_marker_heatmaps_and_myeloid_composition.R
│   ├── 04_slingshot_pseudotime.R
│   ├── 05_cellchat_analysis.R
│   └── 06_final_umap_redraw.R
└── provenance/
    ├── wes/
    └── scrna/
```

The `provenance/` directory contains the selected original scripts or source excerpts from which the cleaned manuscript-facing scripts were prepared. Large data files, RDS objects, FASTQs, BAMs, VCFs and patient material are intentionally not included.

## WES workflow

### 1. Reference preparation

The analysis used GRCh38 (`Homo_sapiens_assembly38.fasta`) and a chromosome-compatible dbSNP known-sites VCF for base-quality score recalibration.

```bash
export REF_DIR=/path/to/reference
bash wes/00_prepare_reference.sh
```

The executed project used a file named `dbsnp_138.hg38.chr.vcf.gz`; the exact upstream dbSNP release was not independently verified from the available VCF header. A separate chromosome-name normalization step was present in the original HPC project because the downloaded dbSNP resource initially used contig names that were incompatible with the reference FASTA.

### 2. Preprocessing and somatic variant calling

```bash
export WORK_DIR=/path/to/WES_GATK
export RAW_DATA=/path/to/WES_fastqs
bash wes/01_preprocess_and_call.sh
```

Per sample, the workflow performs:

1. FastQC
2. paired-end fastp trimming
3. BWA-MEM alignment to GRCh38
4. SAMtools sorting/indexing and flagstat
5. GATK `MarkDuplicates`
6. GATK `BaseRecalibrator` and `ApplyBQSR`
7. GATK `Mutect2` for:
   - `Primary_CIC` vs `Blood`
   - `Recurrent_CIC` vs `Blood`
8. GATK `FilterMutectCalls`
9. `bcftools view -f PASS`
10. MultiQC summary

The fastp settings reproduced from the original analysis are:

- `--detect_adapter_for_pe`
- `--qualified_quality_phred 20`
- `--unqualified_percent_limit 40`
- `--n_base_limit 5`
- `--length_required 50`
- `--correction`

The executed somatic workflow was intentionally simplified: gnomAD, a panel of normals and formal contamination estimation were not supplied to Mutect2, and BQSR used dbSNP as the known-sites resource.

### 3. SnpEff annotation and longitudinal comparison

```bash
export WORK_DIR=/path/to/WES_GATK
bash wes/02_snpeff_annotation_and_compare.sh
bash wes/03_snpeff_parsed_annotation_for_figures.sh
```

SnpEff `GRCh38.99` was used for functional annotation. Primary and recurrent PASS calls are defined as exactly shared only if all four fields are identical:

```text
CHROM + POS + REF + ALT
```

Primary-specific and recurrent-specific calls are set differences between the two PASS call sets. These categories are descriptive and should not be interpreted as formal clonal gain/loss calls.

### 4. WES publication figures

```bash
python wes/generate_wes_publication_figures.py \
    --wes-root /path/to/WES_GATK \
    --output-dir /path/to/figure_output

python wes/generate_exact_shared_vaf_dumbbell.py \
    --wes-root /path/to/WES_GATK \
    --output-dir /path/to/figure_output
```

VAF is recalculated from tumor allele depths as:

```text
ALT_AD / (REF_AD + ALT_AD)
```

No correction is applied for tumor purity, ploidy or local copy-number state. The main figure script contains strict checks for the reviewed publication package; use `--no-strict` only when deliberately applying the code to a different input set.

## scRNA-seq workflow

### 1. Cell Ranger

The recurrent-tumor data were processed with Cell Ranger 8.0.1 and `refdata-gex-GRCh38-2024-A`.

```bash
sbatch scrna/00_cellranger_count.slurm
```

Before submission, set `CELLRANGER_REF` and `FASTQ_DIR`.

### 2. Global Seurat analysis

```bash
Rscript scrna/01_global_seurat_qc_clustering.R \
    /path/to/filtered_feature_bc_matrix \
    results/scRNA/global
```

The manuscript analysis used:

- `CreateSeuratObject(min.cells = 3, min.features = 200)`
- mitochondrial fraction calculated from `^MT-`
- removal of cells with `percent.mt >= 15`
- `LogNormalize`, scale factor 10,000
- 2,000 variable features using `vst`
- PCA followed by neighbors/clustering on PCs 1-20
- global clustering resolution 0.5
- UMAP on PCs 1-20
- positive cluster markers using Wilcoxon testing, `min.pct = 0.25` and `logfc.threshold = 0.25`

No DoubletFinder/Scrublet step or additional upper gene/UMI cutoff is added here because it was not present in the recovered manuscript analysis code.

### 3. Tumor and myeloid reclustering

```bash
Rscript scrna/02_tumor_myeloid_reclustering.R \
    results/scRNA/global/global_major_annotated.rds \
    results/scRNA/reclustered
```

Tumor cells and macrophage/myeloid cells were independently reprocessed with the same general Seurat workflow and reclustered using PCs 1-20 at resolution 0.8.

The tumor compartment was annotated into five states:

- Proliferative tumor
- Mesenchymal/Invasive tumor
- ECM-rich/claudin-high tumor
- Stress-hypoxia epithelial-like tumor
- Neural/Stem-like tumor

The myeloid compartment included resident/homeostatic macrophages, activated TAMs, immunosuppressive/M2-like TAMs, inflammatory myeloid cells, IFN-responsive myeloid cells, neutrophil-like inflammatory myeloid cells, metallothionein/APC-like macrophages and cycling myeloid cells.

### 4. Representative marker heatmaps and myeloid composition

```bash
Rscript scrna/03_marker_heatmaps_and_myeloid_composition.R \
    results/scRNA/reclustered/refined_global.rds \
    results/scRNA/reclustered/Macrophage_Myeloid/Macrophage_Myeloid_reclustered_annotated.rds \
    results/scRNA/markers
```

Representative marker heatmaps use average normalized RNA expression followed by gene-wise z-score scaling. The myeloid composition panel is calculated within the retained annotated myeloid compartment after excluding contaminating/unassigned cells.

### 5. Slingshot pseudotime

```bash
Rscript scrna/04_slingshot_pseudotime.R \
    results/scRNA/reclustered/refined_global.rds \
    results/scRNA/pseudotime
```

The five tumor states are used as Slingshot cluster labels in PCA space. `Neural/Stem-like tumor` is specified as the starting cluster and no terminal state is forced.

The manuscript displays the first inferred lineage. Cells assigned to that lineage are split into pseudotime tertiles (`Early`, `Middle`, `Late`). Stage-associated genes are selected by pairwise Wilcoxon tests with:

- `min.pct = 0.10`
- `logfc.threshold = 0.25`
- adjusted P value `< 0.05`
- positive log2 fold-change against both other stages
- highest average expression in the nominated stage

Up to 40 genes per stage are used for the pseudotime heatmap. Pseudotime is treated as an exploratory representation of transcriptional continuity, not as direct chronological or clonal evolution.

### 6. CellChat

```bash
Rscript scrna/05_cellchat_analysis.R \
    results/scRNA/reclustered/refined_global.rds \
    results/scRNA/cellchat
```

Normalized RNA data are converted from Ensembl identifiers to gene symbols using `org.Hs.eg.db`; unmapped genes are removed and duplicate symbols are summed.

The CellChat workflow uses:

- `CellChatDB.human`
- `subsetData`
- `identifyOverExpressedGenes`
- `identifyOverExpressedInteractions`
- `computeCommunProb`
- `filterCommunication(min.cells = 10)`
- `computeCommunProbPathway`
- `aggregateNet`

Two CellChat runs are retained because both contributed to the manuscript figures:

1. **all refined populations**, used for focused ligand-receptor/pathway ranking;
2. **T/NK and dendritic cells excluded**, used for tumor-to-myeloid and reciprocal myeloid-to-tumor interaction-strength matrices.

The focused manuscript comparison is `Proliferative tumor -> Immunosuppressive / M2-like TAM`.

### 7. Final UMAP rendering

```bash
Rscript scrna/06_final_umap_redraw.R \
    results/scRNA/reclustered/refined_global.rds \
    results/scRNA/reclustered/Tumor_cells/Tumor_reclustered.rds \
    results/scRNA/reclustered/Macrophage_Myeloid/Macrophage_Myeloid_reclustered.rds \
    results/scRNA/final_umaps
```

This script reproduces the final publication display filtering for tumor and myeloid UMAPs. A provenance note about the final myeloid display mapping is provided in `docs/PROVENANCE_AND_LIMITATIONS.md`.

## Reproducibility and interpretation

This repository deliberately preserves several limitations of the executed analysis rather than hiding them:

- WES calling used a simplified Mutect2 configuration.
- WES VAFs are descriptive and are not purity/CNV corrected.
- scRNA-seq cell states were annotated from transcriptional profiles and representative markers; no inferCNV/CopyKAT-based malignancy call is added retrospectively.
- Slingshot is exploratory.
- CellChat reports transcriptome-based communication predictions, not physical receptor engagement or functional signaling.
- This is a single longitudinal case; population-level statistical inference is not appropriate.

See `docs/PROVENANCE_AND_LIMITATIONS.md` for additional details.

## Data

Raw sequencing data and large intermediate objects are not included in this repository. The scripts expect local user-supplied paths to the appropriate FASTQ, reference, VCF or Seurat inputs.

## Provenance

The cleaned scripts were assembled from the original HPC scripts and the manuscript-specific `WHX` section of the recovered R analysis file. Selected source scripts are retained under `provenance/` so that the manuscript-facing code can be traced back to the original analysis.

## Citation

Please cite the associated manuscript when using this analysis code. A full citation will be added after publication.