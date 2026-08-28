# Software versions

This table separates versions that were confirmed by execution/output evidence from versions that appeared only in environment-setup scripts.

## Confirmed from execution evidence

| Software | Version / reference | Evidence |
|---|---|---|
| Cell Ranger | 8.0.1 | Cell Ranger Slurm launcher and run metadata |
| 10x reference | `refdata-gex-GRCh38-2024-A` | Cell Ranger launcher |
| R | 4.3.3 | saved Seurat object metadata extraction |
| Seurat | 5.1.0 | saved Seurat object metadata extraction |
| SeuratObject | 5.0.2 | saved Seurat object metadata extraction |
| FastQC | 0.12.1 | MultiQC software-version output |
| fastp | 1.3.6 | MultiQC software-version output |
| SnpEff | 5.2, build 2023-09-29 | SnpEff execution log |
| SnpEff database | GRCh38.99 | SnpEff execution log |

## Used but exact executed version not recovered

The collected scripts demonstrate use of the following tools, but the final executed version was not reliably recoverable from the available logs:

- BWA-MEM
- SAMtools
- bcftools
- GATK4
- CellChat
- Slingshot
- SingleCellExperiment
- `org.Hs.eg.db`
- ComplexHeatmap

The historical WES environment scripts contain intended/pinned versions (for example BWA 0.7.17, SAMtools/bcftools 1.18, GATK 4.4.0.0, fastp 0.23.4 and SnpEff 5.1), but the observed fastp and SnpEff versions differ from those setup-script values. For that reason, setup-script versions should not be treated as proof of the executed software environment.
