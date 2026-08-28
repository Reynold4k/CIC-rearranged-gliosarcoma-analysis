#!/usr/bin/env Rscript

# Exploratory Slingshot trajectory analysis used for Supplementary Figure S1D-E.
#
# Usage:
#   Rscript 04_slingshot_pseudotime.R \
#       /path/to/refined_global.rds /path/to/output_dir
#
# The manuscript reports the first inferred lineage only. The trajectory is
# inferred in PCA space using the five annotated tumor states, with
# "Neural/Stem-like tumor" specified as the starting cluster and no terminal
# cluster forced.

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_slingshot_pseudotime.R <refined_global.rds> <output_dir>")
}
obj <- readRDS(args[[1]])
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tumor_types <- c(
  "ECM-rich/claudin-high tumor",
  "Proliferative tumor",
  "Mesenchymal/Invasive tumor",
  "Stress-hypoxia epithelial-like tumor",
  "Neural/Stem-like tumor"
)

tumor_cells <- colnames(obj)[obj$celltype_refined %in% tumor_types]
obj_tumor <- subset(obj, cells = tumor_cells)

sce <- as.SingleCellExperiment(obj_tumor)
colData(sce)$slingshot_cluster <- obj_tumor$celltype_refined

sce <- slingshot(
  sce,
  clusterLabels = "slingshot_cluster",
  reducedDim = "PCA",
  start.clus = "Neural/Stem-like tumor"
)

pt_mat <- slingPseudotime(sce)
write.csv(
  as.data.frame(pt_mat) %>% rownames_to_column("cell"),
  file.path(outdir, "slingshot_pseudotime_matrix.csv"),
  row.names = FALSE
)

if (ncol(pt_mat) < 1) {
  stop("Slingshot did not return a lineage.")
}

for (i in seq_len(ncol(pt_mat))) {
  column_name <- paste0("slingshot_pt_lineage", i)
  values <- rep(NA_real_, ncol(obj))
  names(values) <- colnames(obj)
  values[rownames(pt_mat)] <- pt_mat[, i]
  obj[[column_name]] <- values
}
tumor_lineage1 <- rep(NA_real_, ncol(obj_tumor))
names(tumor_lineage1) <- colnames(obj_tumor)
tumor_lineage1[rownames(pt_mat)] <- pt_mat[, 1]
obj_tumor$slingshot_pt_lineage1 <- tumor_lineage1

umap_df <- Embeddings(obj, "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell")
colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")

meta_df <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, celltype_refined, slingshot_pt_lineage1)

plot_df <- umap_df %>%
  left_join(meta_df, by = "cell") %>%
  mutate(is_tumor = celltype_refined %in% tumor_types)

p_pt <- ggplot() +
  geom_point(
    data = plot_df %>% filter(!is_tumor),
    aes(UMAP_1, UMAP_2),
    color = "grey85", size = 0.18, alpha = 0.7
  ) +
  geom_point(
    data = plot_df %>% filter(is_tumor),
    aes(UMAP_1, UMAP_2, color = slingshot_pt_lineage1),
    size = 0.3, alpha = 0.95
  ) +
  scale_color_viridis_c(option = "plasma", na.value = "grey70", begin = 0.05, end = 0.98) +
  coord_equal() +
  theme_classic(base_size = 15) +
  labs(x = "UMAP_1", y = "UMAP_2", color = "Pseudotime") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(linewidth = 0.7, color = "black")
  )

ggsave(file.path(outdir, "S1D_lineage1_pseudotime_global_UMAP.pdf"), p_pt, width = 7.5, height = 6.5)
ggsave(file.path(outdir, "S1D_lineage1_pseudotime_global_UMAP.png"), p_pt, width = 7.5, height = 6.5, dpi = 400)

pt_vec <- obj_tumor$slingshot_pt_lineage1
keep <- !is.na(pt_vec)
obj_lineage <- subset(obj_tumor, cells = colnames(obj_tumor)[keep])
pt_vec <- obj_lineage$slingshot_pt_lineage1

q_cut <- quantile(pt_vec, probs = c(1/3, 2/3), na.rm = TRUE)
obj_lineage$pt_stage <- factor(
  dplyr::case_when(
    pt_vec <= q_cut[1] ~ "Early",
    pt_vec > q_cut[1] & pt_vec <= q_cut[2] ~ "Middle",
    pt_vec > q_cut[2] ~ "Late"
  ),
  levels = c("Early", "Middle", "Late")
)

write.csv(
  obj_lineage@meta.data %>%
    rownames_to_column("cell") %>%
    dplyr::select(cell, celltype_refined, slingshot_pt_lineage1, pt_stage),
  file.path(outdir, "S1E_lineage1_cells_with_stage.csv"),
  row.names = FALSE
)

DefaultAssay(obj_lineage) <- "RNA"
Idents(obj_lineage) <- obj_lineage$pt_stage

find_pair <- function(a, b) {
  FindMarkers(
    obj_lineage,
    ident.1 = a,
    ident.2 = b,
    assay = "RNA",
    logfc.threshold = 0.25,
    min.pct = 0.10,
    test.use = "wilcox"
  ) %>% rownames_to_column("gene")
}

deg_E_M <- find_pair("Early", "Middle")
deg_E_L <- find_pair("Early", "Late")
deg_M_E <- find_pair("Middle", "Early")
deg_M_L <- find_pair("Middle", "Late")
deg_L_E <- find_pair("Late", "Early")
deg_L_M <- find_pair("Late", "Middle")

write.csv(deg_E_M, file.path(outdir, "Early_vs_Middle_DEG.csv"), row.names = FALSE)
write.csv(deg_E_L, file.path(outdir, "Early_vs_Late_DEG.csv"), row.names = FALSE)
write.csv(deg_M_E, file.path(outdir, "Middle_vs_Early_DEG.csv"), row.names = FALSE)
write.csv(deg_M_L, file.path(outdir, "Middle_vs_Late_DEG.csv"), row.names = FALSE)
write.csv(deg_L_E, file.path(outdir, "Late_vs_Early_DEG.csv"), row.names = FALSE)
write.csv(deg_L_M, file.path(outdir, "Late_vs_Middle_DEG.csv"), row.names = FALSE)

early_genes <- intersect(
  deg_E_M %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene),
  deg_E_L %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene)
)
middle_genes <- intersect(
  deg_M_E %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene),
  deg_M_L %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene)
)
late_genes <- intersect(
  deg_L_E %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene),
  deg_L_M %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% pull(gene)
)

avg_expr <- AverageExpression(
  obj_lineage,
  group.by = "pt_stage",
  assays = "RNA",
  layer = "data",
  return.seurat = FALSE
)$RNA
avg_expr_df <- as.data.frame(avg_expr) %>% rownames_to_column("gene")

early_genes <- avg_expr_df %>%
  filter(gene %in% early_genes, Early > Middle, Early > Late) %>%
  arrange(desc(Early - pmax(Middle, Late))) %>%
  pull(gene)
middle_genes <- avg_expr_df %>%
  filter(gene %in% middle_genes, Middle > Early, Middle > Late) %>%
  arrange(desc(Middle - pmax(Early, Late))) %>%
  pull(gene)
late_genes <- avg_expr_df %>%
  filter(gene %in% late_genes, Late > Early, Late > Middle) %>%
  arrange(desc(Late - pmax(Early, Middle))) %>%
  pull(gene)

write.csv(data.frame(gene = early_genes), file.path(outdir, "Early_specific_genes.csv"), row.names = FALSE)
write.csv(data.frame(gene = middle_genes), file.path(outdir, "Middle_specific_genes.csv"), row.names = FALSE)
write.csv(data.frame(gene = late_genes), file.path(outdir, "Late_specific_genes.csv"), row.names = FALSE)

top_early <- head(early_genes, 40)
top_middle <- head(middle_genes, 40)
top_late <- head(late_genes, 40)
heat_genes <- unique(c(top_early, top_middle, top_late))

cells_use <- rownames(obj_lineage@meta.data)
cells_use <- cells_use[order(obj_lineage@meta.data[cells_use, "slingshot_pt_lineage1"])]

expr_mat <- GetAssayData(obj_lineage, assay = "RNA", layer = "data")[heat_genes, cells_use, drop = FALSE]
expr_mat <- as.matrix(expr_mat)
expr_scaled <- t(scale(t(expr_mat)))
expr_scaled[is.na(expr_scaled)] <- 0
expr_scaled[expr_scaled > 2.5] <- 2.5
expr_scaled[expr_scaled < -2.5] <- -2.5

gene_group <- c(
  rep("Early", length(top_early)),
  rep("Middle", length(top_middle)),
  rep("Late", length(top_late))
)
names(gene_group) <- c(top_early, top_middle, top_late)
gene_group <- factor(gene_group[rownames(expr_scaled)], levels = c("Early", "Middle", "Late"))

cell_stage <- obj_lineage$pt_stage[match(colnames(expr_scaled), colnames(obj_lineage))]
pt_range <- range(obj_lineage$slingshot_pt_lineage1, na.rm = TRUE)

ha <- HeatmapAnnotation(
  Stage = cell_stage,
  Pseudotime = obj_lineage$slingshot_pt_lineage1[match(colnames(expr_scaled), colnames(obj_lineage))],
  col = list(
    Stage = c(Early = "#8FB7A3", Middle = "#C2A27A", Late = "#A99BC2"),
    Pseudotime = colorRamp2(
      c(pt_range[1], mean(pt_range), pt_range[2]),
      c("#440154", "#FDE725", "#B12A90")
    )
  )
)

ht <- Heatmap(
  expr_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2.5, 0, 2.5), c("#3B4CC0", "white", "#B40426")),
  top_annotation = ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  show_column_names = FALSE,
  row_split = gene_group,
  use_raster = TRUE,
  raster_quality = 3,
  border = FALSE,
  column_title = "Lineage 1 cells ordered by pseudotime"
)

pdf(file.path(outdir, "S1E_lineage1_stage_specific_heatmap.pdf"), width = 12, height = 10, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(file.path(outdir, "S1E_lineage1_stage_specific_heatmap.png"), width = 12, height = 10, units = "in", res = 400)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

saveRDS(obj, file.path(outdir, "refined_global_with_slingshot.rds"))
