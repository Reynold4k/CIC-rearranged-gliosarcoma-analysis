#!/usr/bin/env Rscript

# Global Seurat preprocessing, QC, clustering, marker detection, and major
# cell-type annotation used for the recurrent tumor scRNA-seq analysis.
#
# Usage:
#   Rscript 01_global_seurat_qc_clustering.R \
#       /path/to/filtered_feature_bc_matrix /path/to/output_dir
#
# Analysis parameters are reproduced from the original manuscript workflow.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 01_global_seurat_qc_clustering.R <filtered_matrix_dir> <output_dir>")
}
data_dir <- normalizePath(args[[1]], mustWork = TRUE)
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts <- Read10X(data.dir = data_dir)

obj <- CreateSeuratObject(
  counts = counts,
  project = "recurrent_tumor_cellranger",
  min.cells = 3,
  min.features = 200
)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

pdf(file.path(outdir, "QC_violin_before_filtering.pdf"), width = 10, height = 4)
print(VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
))
dev.off()

before_n <- ncol(obj)
obj <- subset(obj, subset = percent.mt < 15)
after_n <- ncol(obj)
writeLines(
  c(
    paste("Before filtering:", before_n),
    paste("After filtering:", after_n),
    paste("Removed:", before_n - after_n)
  ),
  file.path(outdir, "QC_cell_counts.txt")
)

obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2000
)
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

pdf(file.path(outdir, "ElbowPlot.pdf"), width = 6, height = 5)
print(ElbowPlot(obj))
dev.off()

obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:20)

markers_all <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
top10_markers <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

write.csv(markers_all, file.path(outdir, "all_cluster_markers.csv"), row.names = FALSE)
write.csv(top10_markers, file.path(outdir, "top10_cluster_markers.csv"), row.names = FALSE)

# Manual major-lineage annotation used in the source analysis.
cluster2major <- c(
  "0"  = "Tumor cells",
  "1"  = "Macrophage / Myeloid",
  "2"  = "Macrophage / Myeloid",
  "3"  = "Tumor cells",
  "4"  = "Tumor cells",
  "5"  = "T / NK cells",
  "6"  = "Macrophage / Myeloid",
  "7"  = "Tumor cells",
  "8"  = "Macrophage / Myeloid",
  "9"  = "Platelet / Megakaryocyte",
  "10" = "Fibroblast (CAF)",
  "11" = "Macrophage / Myeloid",
  "12" = "Endothelial cells",
  "13" = "Neuron",
  "14" = "Dendritic cells",
  "15" = "Dendritic cells",
  "16" = "Tumor cells",
  "17" = "Endothelial cells"
)
obj$celltype_major <- unname(cluster2major[as.character(Idents(obj))])

p_major <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "celltype_major",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  pt.size = 0.6,
  raster = FALSE
) +
  ggtitle("Major cell types") +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(outdir, "UMAP_celltype_major.pdf"), p_major, width = 8, height = 7)
ggsave(file.path(outdir, "UMAP_celltype_major.png"), p_major, width = 8, height = 7, dpi = 600)

saveRDS(obj, file.path(outdir, "global_major_annotated.rds"))
