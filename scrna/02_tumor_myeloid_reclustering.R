#!/usr/bin/env Rscript

# Tumor and macrophage/myeloid reclustering, subtype annotation, and construction
# of the refined global annotation used by pseudotime and CellChat.
#
# Usage:
#   Rscript 02_tumor_myeloid_reclustering.R \
#       /path/to/global_major_annotated.rds /path/to/output_dir

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 02_tumor_myeloid_reclustering.R <global_rds> <output_dir>")
}
obj <- readRDS(args[[1]])
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
tumor_outdir <- file.path(outdir, "Tumor_cells")
myeloid_outdir <- file.path(outdir, "Macrophage_Myeloid")
dir.create(tumor_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(myeloid_outdir, recursive = TRUE, showWarnings = FALSE)

recluster_and_save <- function(seu, outdir, prefix, dims_use = 1:20, resolution = 0.8) {
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)
  seu <- FindNeighbors(seu, dims = dims_use, verbose = FALSE)
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  seu <- RunUMAP(seu, dims = dims_use, verbose = FALSE)

  markers_all <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  top10_markers <- markers_all %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
    ungroup()

  write.csv(markers_all, file.path(outdir, paste0(prefix, "_all_cluster_markers.csv")), row.names = FALSE)
  write.csv(top10_markers, file.path(outdir, paste0(prefix, "_top10_cluster_markers.csv")), row.names = FALSE)

  p <- DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE) +
    ggtitle(paste0(prefix, " reclustering"))
  ggsave(file.path(outdir, paste0(prefix, "_UMAP_recluster.pdf")), p, width = 8, height = 6)
  ggsave(file.path(outdir, paste0(prefix, "_UMAP_recluster.png")), p, width = 8, height = 6, dpi = 300)

  saveRDS(seu, file.path(outdir, paste0(prefix, "_reclustered.rds")))
  seu
}

tumor_obj <- subset(obj, subset = celltype_major == "Tumor cells")
tumor_obj <- recluster_and_save(tumor_obj, tumor_outdir, "Tumor", 1:20, 0.8)

macro_obj <- subset(obj, subset = celltype_major == "Macrophage / Myeloid")
macro_obj <- recluster_and_save(macro_obj, myeloid_outdir, "Macrophage_Myeloid", 1:20, 0.8)

# Tumor-state labels used by the downstream manuscript analysis.
tumor_cluster2type <- c(
  "0" = "ECM-rich/claudin-high tumor",
  "1" = "Proliferative tumor",
  "2" = "Mesenchymal/Invasive tumor",
  "4" = "Neural/Stem-like tumor",
  "5" = "Proliferative tumor",
  "6" = "Stress-hypoxia epithelial-like tumor",
  "7" = "Neural/Stem-like tumor",
  "8" = "Proliferative tumor"
)
tumor_obj$tumor_subtype <- unname(tumor_cluster2type[as.character(Idents(tumor_obj))])
tumor_obj$tumor_subtype[is.na(tumor_obj$tumor_subtype)] <- "Exclude_from_tumor"

# Myeloid-state labels used in the original downstream analysis object.
# The final publication UMAP was later redrawn with a refined display mapping;
# see 06_final_umap_redraw.R and docs/PROVENANCE_AND_LIMITATIONS.md.
myeloid_cluster2type_analysis <- c(
  "0"  = "Resident/homeostatic macrophage",
  "1"  = "Activated TAM",
  "2"  = "Immunosuppressive / M2-like TAM",
  "3"  = "Immunosuppressive / M2-like TAM",
  "4"  = "Cycling myeloid",
  "5"  = "Inflammatory myeloid",
  "6"  = "Metallothionein / APC-like macrophage",
  "7"  = "Cycling myeloid",
  "8"  = "Neutrophil-like inflammatory myeloid",
  "9"  = "IFN-responsive myeloid",
  "10" = "T/NK contamination",
  "11" = "Cycling myeloid"
)
macro_obj$myeloid_subtype_final <- unname(
  myeloid_cluster2type_analysis[as.character(Idents(macro_obj))]
)
macro_obj$myeloid_subtype_final[is.na(macro_obj$myeloid_subtype_final)] <- "Unassigned"

# Construct the refined global annotation used in the original pseudotime and
# CellChat workflow.
obj <- subset(obj, subset = celltype_major != "Platelet / Megakaryocyte")
obj$celltype_refined <- obj$celltype_major
obj$celltype_refined[colnames(macro_obj)] <- macro_obj$myeloid_subtype_final
obj$celltype_refined[colnames(tumor_obj)] <- tumor_obj$tumor_subtype
obj$celltype_refined[obj$celltype_refined == "Exclude_from_tumor"] <- "myeloid/macrophage/TAM"
obj$celltype_refined <- dplyr::recode(
  obj$celltype_refined,
  "Macrophage / Myeloid" = "myeloid/macrophage/TAM"
)

saveRDS(tumor_obj, file.path(tumor_outdir, "Tumor_reclustered_annotated.rds"))
saveRDS(macro_obj, file.path(myeloid_outdir, "Macrophage_Myeloid_reclustered_annotated.rds"))
saveRDS(obj, file.path(outdir, "refined_global.rds"))
