# Methods evidence summary

## WES

### Preprocessing
The recovered WES pipeline applies FastQC, paired-end fastp trimming, BWA-MEM alignment to GRCh38, SAMtools sorting/indexing, GATK MarkDuplicates and GATK base-quality score recalibration.

### Somatic calling
Primary and recurrent tumors are each called against the same blood normal with GATK Mutect2. FilterMutectCalls is applied and only PASS records are retained.

### Annotation and longitudinal comparison
SnpEff GRCh38.99 is used for functional annotation. Shared variants are defined by exact `CHROM + POS + REF + ALT` identity. VAF in the figure scripts is derived from tumor AD values as `ALT/(REF+ALT)` and is not purity/CNV corrected.

## scRNA-seq

### Cell Ranger
The recurrent tumor dataset was processed with Cell Ranger 8.0.1 using `refdata-gex-GRCh38-2024-A`, 24 local cores and 128 GB local memory.

### Seurat QC and global clustering
Recovered code shows `min.cells=3`, `min.features=200`, `percent.mt < 15`, LogNormalize with scale factor 10,000, 2,000 VST variable genes, PCA, PCs 1-20 for neighbors/clustering, global resolution 0.5 and Wilcoxon positive-marker detection (`min.pct=0.25`, `logfc.threshold=0.25`).

### Tumor/myeloid reclustering
Both compartments are independently normalized, variable-feature selected, scaled, subjected to PCA, graph clustering and UMAP using PCs 1-20 and resolution 0.8.

### Marker heatmaps
Representative marker genes are summarized by average normalized RNA expression per annotated population and scaled by gene (row-wise z-score).

### Slingshot
Five tumor states are used as cluster labels in PCA space. `Neural/Stem-like tumor` is the starting cluster. The manuscript uses the first inferred lineage. Pseudotime is divided into tertiles and stage-associated genes are identified by pairwise Wilcoxon tests with the thresholds implemented in `04_slingshot_pseudotime.R`.

### CellChat
Normalized RNA expression is converted from Ensembl IDs to gene symbols; duplicate symbols are combined. `CellChatDB.human` and `filterCommunication(min.cells=10)` are used. Both the all-cell and no-T/NK/no-dendritic runs that contributed to figures are retained.
