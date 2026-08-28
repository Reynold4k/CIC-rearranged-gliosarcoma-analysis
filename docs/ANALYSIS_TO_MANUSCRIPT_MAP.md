# Analysis-to-manuscript map

This document maps repository scripts to the analyses described in the manuscript and its figure legends.

| Manuscript element | Analysis | Primary code |
|---|---|---|
| Figure 1 WES panel(s) | exact shared PASS variants and AD-based VAF visualization | `wes/generate_exact_shared_vaf_dumbbell.py`, `wes/generate_wes_publication_figures.py` |
| Supplementary Fig. S3A | WES Q30/mapping/duplication summaries | `wes/01_preprocess_and_call.sh` plus MultiQC |
| Supplementary Fig. S3B-C | PASS counts and exact primary/shared/recurrent membership | `wes/02_snpeff_annotation_and_compare.sh`, `wes/generate_wes_publication_figures.py` |
| Supplementary Fig. S3D | primary vs recurrent AD-based VAFs for exact shared variants | `wes/generate_wes_publication_figures.py` |
| Figure 2A | global single-cell UMAP / major populations | `scrna/01_global_seurat_qc_clustering.R` |
| Figure 2B | tumor-cell state UMAP | `scrna/02_tumor_myeloid_reclustering.R`, `scrna/06_final_umap_redraw.R` |
| Figure 2C | macrophage/myeloid state UMAP | `scrna/02_tumor_myeloid_reclustering.R`, `scrna/06_final_umap_redraw.R` |
| Figure 2D | predicted tumor-to-myeloid interaction-strength matrix | `scrna/05_cellchat_analysis.R` (round 2) |
| Figure 2E | focused proliferative-tumor -> M2-like TAM ligand-receptor ranking | `scrna/05_cellchat_analysis.R` (round 1) |
| Figure 2F | focused proliferative-tumor -> M2-like TAM pathway ranking | `scrna/05_cellchat_analysis.R` (round 1) |
| Supplementary Fig. S1A | major-cell marker heatmap | `scrna/03_marker_heatmaps_and_myeloid_composition.R` |
| Supplementary Fig. S1B | tumor-state marker heatmap | `scrna/03_marker_heatmaps_and_myeloid_composition.R` |
| Supplementary Fig. S1C | refined global annotation UMAP | `scrna/06_final_umap_redraw.R` |
| Supplementary Fig. S1D | lineage-1 pseudotime projected on global UMAP | `scrna/04_slingshot_pseudotime.R` |
| Supplementary Fig. S1E | pseudotime-associated Early/Middle/Late heatmap | `scrna/04_slingshot_pseudotime.R` |
| Supplementary Fig. S2A | relative myeloid-state composition | `scrna/03_marker_heatmaps_and_myeloid_composition.R` |
| Supplementary Fig. S2B | myeloid-state marker heatmap | `scrna/03_marker_heatmaps_and_myeloid_composition.R` |
| Supplementary Fig. S2C | broader predicted pathway relationships | `scrna/05_cellchat_analysis.R` (round 1 priority-pair pathway overview) |
| Supplementary Fig. S2D | predicted myeloid-to-tumor interaction strengths | `scrna/05_cellchat_analysis.R` (round 2) |

Multiplex immunofluorescence in Figure 2G is not a computational analysis and is not represented by an analysis script in this repository.
