#!/usr/bin/env Rscript

# Final publication UMAP rendering for tumor and myeloid subtypes, plus the
# refined global annotation shown in Supplementary Figure S1C.
#
# Usage:
#   Rscript 06_final_umap_redraw.R \
#       /path/to/refined_global.rds \
#       /path/to/Tumor_reclustered.rds \
#       /path/to/Macrophage_Myeloid_reclustered.rds \
#       /path/to/output_dir
#
# IMPORTANT: the final myeloid display mapping in the original redraw script
# differed from the earlier exploratory mapping used to construct the
# downstream refined analysis object. This script reproduces the final display
# mapping exactly. See docs/PROVENANCE_AND_LIMITATIONS.md.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage: Rscript 06_final_umap_redraw.R ",
    "<refined_global.rds> <tumor_reclustered.rds> <myeloid_reclustered.rds> <output_dir>"
  )
}

refined_obj <- readRDS(args[[1]])
tumor_obj <- readRDS(args[[2]])
myeloid_obj <- readRDS(args[[3]])
outdir <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

save_pair <- function(plot, stem) {
  ggsave(file.path(outdir, paste0(stem, ".png")), plot, width = 9, height = 7, dpi = 600, bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".pdf")), plot, width = 9, height = 7, device = cairo_pdf, bg = "white")
}

# Supplementary Figure S1C: refined global annotation.
p_global <- DimPlot(
  refined_obj,
  reduction = "umap",
  group.by = "celltype_refined",
  label = TRUE,
  repel = TRUE,
  label.size = 5.5,
  pt.size = 0.6,
  raster = FALSE
) +
  labs(title = "Refined cell type annotation", color = NULL) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    aspect.ratio = 1
  )
save_pair(p_global, "S1C_refined_global_UMAP")

# Figure 2B: final tumor-state UMAP. Clusters 3, 9 and 10 were excluded from
# the tumor-specific display in the original final plotting script.
tumor_map <- c(
  "0" = "ECM-rich/claudin-high tumor",
  "1" = "Proliferative tumor",
  "2" = "Mesenchymal/Invasive tumor",
  "3" = "Exclude_from_tumor",
  "4" = "Neural/Stem-like tumor",
  "5" = "Proliferative tumor",
  "6" = "Stress-hypoxia epithelial-like tumor",
  "7" = "Neural/Stem-like tumor",
  "8" = "Proliferative tumor",
  "9" = "Exclude_from_tumor",
  "10" = "Exclude_from_tumor"
)
tumor_subtype <- unname(tumor_map[as.character(tumor_obj$seurat_clusters)])
keep_tumor <- tumor_subtype != "Exclude_from_tumor"
tumor_obj$tumor_subtype_plot <- factor(
  ifelse(keep_tumor, tumor_subtype, NA_character_),
  levels = c(
    "ECM-rich/claudin-high tumor",
    "Mesenchymal/Invasive tumor",
    "Neural/Stem-like tumor",
    "Proliferative tumor",
    "Stress-hypoxia epithelial-like tumor"
  )
)
p_tumor <- DimPlot(
  tumor_obj,
  reduction = "umap",
  cells = colnames(tumor_obj)[keep_tumor],
  group.by = "tumor_subtype_plot",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.5
) +
  labs(title = "Tumor subtype annotation", color = NULL) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_pair(p_tumor, "Fig2B_tumor_subtype_UMAP")

# Figure 2C: final myeloid display mapping from redraw_filtered_umaps.R.
myeloid_map <- c(
  "0" = "Resident/homeostatic macrophage",
  "1" = "Activated TAM",
  "2" = "IFN-responsive myeloid",
  "3" = "Immunosuppressive / M2-like TAM",
  "4" = "Cycling myeloid",
  "5" = "Inflammatory myeloid",
  "6" = "Metallothionein / APC-like macrophage",
  "7" = "Cycling myeloid",
  "8" = "Neutrophil-like inflammatory myeloid",
  "9" = "IFN-responsive myeloid",
  "10" = "T/NK contamination",
  "11" = "T/NK contamination"
)
myeloid_subtype <- unname(myeloid_map[as.character(myeloid_obj$seurat_clusters)])
keep_myeloid <- myeloid_subtype != "T/NK contamination"
myeloid_obj$myeloid_subtype_plot <- factor(
  ifelse(keep_myeloid, myeloid_subtype, NA_character_),
  levels = c(
    "Activated TAM",
    "Cycling myeloid",
    "IFN-responsive myeloid",
    "Immunosuppressive / M2-like TAM",
    "Inflammatory myeloid",
    "Metallothionein / APC-like macrophage",
    "Neutrophil-like inflammatory myeloid",
    "Resident/homeostatic macrophage"
  )
)
p_myeloid <- DimPlot(
  myeloid_obj,
  reduction = "umap",
  cells = colnames(myeloid_obj)[keep_myeloid],
  group.by = "myeloid_subtype_plot",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.5
) +
  labs(title = "Macrophage / Myeloid subtype annotation", color = NULL) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_pair(p_myeloid, "Fig2C_myeloid_subtype_UMAP")
