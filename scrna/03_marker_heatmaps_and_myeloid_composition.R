#!/usr/bin/env Rscript

# Representative marker heatmaps and the myeloid-composition plot used in
# Supplementary Figures S1 and S2.
#
# Usage:
#   Rscript 03_marker_heatmaps_and_myeloid_composition.R \
#       /path/to/refined_global.rds \
#       /path/to/Macrophage_Myeloid_reclustered_annotated.rds \
#       /path/to/output_dir

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript 03_marker_heatmaps_and_myeloid_composition.R ",
    "<refined_global.rds> <myeloid_rds> <output_dir>"
  )
}

obj <- readRDS(args[[1]])
macro_obj <- readRDS(args[[2]])
outdir <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
DefaultAssay(obj) <- "RNA"

major_marker_list <- list(
  "Tumor cells" = c("HMGA2","PDGFRA","COL6A1","CLDN7"),
  "myeloid/macrophage/TAM" = c("LYZ","FCER1G","C1QC","APOE"),
  "T / NK cells" = c("CD3D","NKG7","GNLY","PRF1"),
  "Dendritic cells" = c("FCER1A","HLA-DRA","CLEC9A","XCR1"),
  "Fibroblast (CAF)" = c("COL1A1","DCN","LUM","POSTN"),
  "Endothelial cells" = c("VWF","EMCN","CLDN5","ROBO4"),
  "Neuron" = c("NLGN1","PRIMA1","SNAP25","RBFOX3")
)

tumor_marker_list <- list(
  "Proliferative tumor" = c("MKI67","TOP2A","UBE2C","CDC20","CCNB1"),
  "Mesenchymal/Invasive tumor" = c("IGFBP5","FRAS1","MEOX2","FN1","COL1A1"),
  "ECM-rich/claudin-high tumor" = c("CLDN1","CLDN7","DCN","BMP3","CA12"),
  "Stress-hypoxia epithelial-like tumor" = c("KRT17","KRT18","NUPR1","HIGD1A","LDHA"),
  "Neural/Stem-like tumor" = c("HMGA2","PDGFRA","NES","SOX2","GPC5")
)

myeloid_marker_list <- list(
  "Resident/homeostatic macrophage" = c("SELENOP","IGF1","TMEM176B","SLC40A1"),
  "Immunosuppressive / M2-like TAM" = c("CCL18","APOE","APOC1","TREM2","CXCL16"),
  "Activated TAM" = c("OLR1","EREG","ITGA5","ANPEP","THBS1"),
  "Metallothionein / APC-like macrophage" = c("MT1X","MT1E","ABI3","CAMK1","HLA-DRA"),
  "Inflammatory myeloid" = c("FCN1","VCAN","IL1B","CXCL8","LYZ"),
  "IFN-responsive myeloid" = c("CXCL10","CXCL11","RSAD2","IFIT1","IFIT3","ISG15"),
  "Neutrophil-like inflammatory myeloid" = c("FCGR3B","S100A12","S100A8","S100A9","IL1R2"),
  "Cycling myeloid" = c("MKI67","TOP2A","UBE2C","CDC20","PLK1","CCNB1")
)

celltype_order <- c(
  "Tumor cells",
  names(tumor_marker_list),
  "myeloid/macrophage/TAM",
  names(myeloid_marker_list),
  "T / NK cells",
  "Dendritic cells",
  "Fibroblast (CAF)",
  "Endothelial cells",
  "Neuron"
)

obj_clean <- subset(obj, subset = !is.na(celltype_refined))
obj_clean$celltype_refined <- droplevels(factor(obj_clean$celltype_refined))
order_use <- celltype_order[celltype_order %in% unique(as.character(obj_clean$celltype_refined))]
obj_clean$celltype_refined <- factor(obj_clean$celltype_refined, levels = order_use)

make_marker_df <- function(marker_list, class_name) {
  do.call(rbind, lapply(names(marker_list), function(x) {
    data.frame(
      module = x,
      gene = marker_list[[x]],
      class = class_name,
      stringsAsFactors = FALSE
    )
  }))
}

marker_df <- bind_rows(
  make_marker_df(major_marker_list, "Major"),
  make_marker_df(tumor_marker_list, "Tumor"),
  make_marker_df(myeloid_marker_list, "Myeloid")
) %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(gene %in% rownames(obj_clean))

avg_expr <- AverageExpression(
  obj_clean,
  features = marker_df$gene,
  group.by = "celltype_refined",
  assays = "RNA",
  layer = "data",
  verbose = FALSE
)$RNA

avg_expr <- avg_expr[marker_df$gene, order_use, drop = FALSE]
mat_scaled <- t(scale(t(avg_expr)))
mat_scaled[is.na(mat_scaled)] <- 0

row_anno <- marker_df[, c("gene", "module", "class")]
rownames(row_anno) <- row_anno$gene
row_anno$gene <- NULL
row_anno <- row_anno[rownames(mat_scaled), , drop = FALSE]

# Deterministic display colors. These affect annotation bars only, not values.
module_names <- unique(as.character(row_anno$module))
display_cols <- grDevices::hcl.colors(length(module_names), "Dark 3")
module_col_map <- setNames(display_cols, module_names)

plot_module_heatmap <- function(modules_keep, out_pdf, title, width, height, row_fontsize = 9) {
  keep <- rownames(row_anno)[as.character(row_anno$module) %in% modules_keep]
  mat_sub <- mat_scaled[keep, , drop = FALSE]
  anno_sub <- row_anno[keep, , drop = FALSE]
  anno_sub$module <- factor(as.character(anno_sub$module), levels = modules_keep)
  ord <- order(anno_sub$module)
  anno_sub <- anno_sub[ord, , drop = FALSE]
  mat_sub <- mat_sub[rownames(anno_sub), , drop = FALSE]

  row_ha <- rowAnnotation(
    Module = anno_sub$module,
    col = list(Module = module_col_map[unique(as.character(anno_sub$module))]),
    show_annotation_name = FALSE,
    width = unit(4, "mm")
  )

  ht <- Heatmap(
    mat_sub,
    name = "z-score",
    col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = row_fontsize),
    column_names_gp = gpar(fontsize = 11, fontface = "bold"),
    row_split = anno_sub$module,
    left_annotation = row_ha,
    row_title = NULL,
    row_gap = unit(4, "mm"),
    column_title = title
  )

  pdf(out_pdf, width = width, height = height)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

plot_module_heatmap(
  names(major_marker_list),
  file.path(outdir, "S1A_major_cell_marker_heatmap.pdf"),
  "Major cell type marker heatmap",
  16, 12, 10
)
plot_module_heatmap(
  names(tumor_marker_list),
  file.path(outdir, "S1B_tumor_state_marker_heatmap.pdf"),
  "Tumor subtype marker heatmap",
  16, 14, 10
)
plot_module_heatmap(
  names(myeloid_marker_list),
  file.path(outdir, "S2B_myeloid_state_marker_heatmap.pdf"),
  "Myeloid subtype marker heatmap",
  16, 18, 9
)

# Supplementary Figure S2A: composition within retained annotated myeloid cells.
myeloid_order <- names(myeloid_marker_list)
x <- as.character(macro_obj$myeloid_subtype_final)
x <- x[!is.na(x) & x != "" & !x %in% c("T/NK contamination", "Unassigned")]
df <- as.data.frame(table(x), stringsAsFactors = FALSE)
colnames(df) <- c("celltype", "n")
df$proportion <- df$n / sum(df$n)
df$percent <- 100 * df$proportion
df$celltype <- factor(df$celltype, levels = myeloid_order[myeloid_order %in% df$celltype])

fill_map <- setNames(
  grDevices::hcl.colors(length(levels(droplevels(df$celltype))), "Dark 3"),
  levels(droplevels(df$celltype))
)

p <- ggplot(df, aes(x = 1, y = proportion, fill = celltype)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent)),
    position = position_stack(vjust = 0.5),
    size = 3.2,
    fontface = "bold"
  ) +
  scale_fill_manual(values = fill_map, name = "Myeloid subtype") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Myeloid subtypes", x = NULL, y = "Proportion") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

ggsave(file.path(outdir, "S2A_myeloid_subtype_composition.pdf"), p, width = 7, height = 10)
ggsave(file.path(outdir, "S2A_myeloid_subtype_composition.png"), p, width = 7, height = 10, dpi = 300)
write.csv(df, file.path(outdir, "S2A_myeloid_subtype_composition.csv"), row.names = FALSE)
