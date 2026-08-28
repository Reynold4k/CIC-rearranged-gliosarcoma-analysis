#!/usr/bin/env Rscript

# CellChat analysis used for Figure 2D-F and Supplementary Figure S2C-D.
#
# Usage:
#   Rscript 05_cellchat_analysis.R /path/to/refined_global.rds /path/to/output_dir
#
# Two rounds are reproduced from the original workflow:
#   1) all refined cell types, used for pathway/LR ranking;
#   2) T/NK and dendritic cells removed, used for tumor<->myeloid interaction
#      strength matrices.
#
# Interactions involving cell groups with fewer than 10 cells are filtered.

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(forcats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 05_cellchat_analysis.R <refined_global.rds> <output_dir>")
}
obj <- readRDS(args[[1]])
outdir <- args[[2]]
round1_dir <- file.path(outdir, "round1_all")
round2_dir <- file.path(outdir, "round2_no_TNK_DC")
dir.create(round1_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(round2_dir, recursive = TRUE, showWarnings = FALSE)

prepare_cellchat_input <- function(seu) {
  DefaultAssay(seu) <- "RNA"
  data_input <- GetAssayData(seu, assay = "RNA", layer = "data")

  original_ids <- rownames(data_input)
  clean_ids <- sub("\\..*$", "", original_ids)
  is_ensg <- grepl("^ENSG[0-9]+$", clean_ids)

  mapped <- rep(NA_character_, length(clean_ids))
  mapped[is_ensg] <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = clean_ids[is_ensg],
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  mapped[!is_ensg] <- clean_ids[!is_ensg]

  keep <- !is.na(mapped) & mapped != ""
  mat <- data_input[keep, ]
  rownames(mat) <- mapped[keep]

  mat <- rowsum(as.matrix(mat), group = rownames(mat), reorder = FALSE)
  as(mat, "dgCMatrix")
}

run_cellchat <- function(seu) {
  mat <- prepare_cellchat_input(seu)
  meta <- data.frame(
    celltype = seu$celltype_refined,
    row.names = colnames(seu)
  )
  meta <- meta[!is.na(meta$celltype) & meta$celltype != "", , drop = FALSE]
  mat <- mat[, rownames(meta), drop = FALSE]

  cc <- createCellChat(object = mat, meta = meta, group.by = "celltype")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# Round 1: all refined cell types
cellchat <- run_cellchat(obj)
saveRDS(cellchat, file.path(round1_dir, "cellchat_all_refined.rds"))

lr_all <- subsetCommunication(cellchat)
pathway_all <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(lr_all, file.path(round1_dir, "CellChat_all_LR_communications.csv"), row.names = FALSE)
write.csv(pathway_all, file.path(round1_dir, "CellChat_pathway_level_communications.csv"), row.names = FALSE)
write.csv(cellchat@net$count, file.path(round1_dir, "CellChat_interaction_count_matrix.csv"))
write.csv(cellchat@net$weight, file.path(round1_dir, "CellChat_interaction_weight_matrix.csv"))

priority_pairs <- data.frame(
  source = c(
    "Proliferative tumor",
    "Proliferative tumor",
    "Mesenchymal/Invasive tumor",
    "Fibroblast (CAF)"
  ),
  target = c(
    "Immunosuppressive / M2-like TAM",
    "Endothelial cells",
    "Fibroblast (CAF)",
    "Proliferative tumor"
  ),
  pair_label = c(
    "Prolif. tumor -> M2-TAM",
    "Prolif. tumor -> Endothelial",
    "Mesench. tumor -> CAF",
    "CAF -> Prolif. tumor"
  ),
  stringsAsFactors = FALSE
)

bubble_df <- pathway_all %>%
  inner_join(priority_pairs, by = c("source", "target")) %>%
  group_by(source, target, pair_label) %>%
  arrange(desc(prob), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

if (nrow(bubble_df) > 0) {
  pathway_order <- bubble_df %>%
    group_by(pathway_name) %>%
    summarise(total_prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(total_prob) %>%
    pull(pathway_name)

  bubble_df$pair_label <- factor(bubble_df$pair_label, levels = priority_pairs$pair_label)
  bubble_df$pathway_name <- factor(bubble_df$pathway_name, levels = pathway_order)

  p_bubble <- ggplot(
    bubble_df,
    aes(x = pair_label, y = pathway_name, size = prob, color = prob)
  ) +
    geom_point(alpha = 0.9) +
    scale_color_viridis_c(option = "D", name = "Probability") +
    scale_size_continuous(name = "Probability") +
    labs(x = NULL, y = "Signaling pathway") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(round1_dir, "S2C_priority_pair_pathways.pdf"), p_bubble, width = 10, height = 8)
  ggsave(file.path(round1_dir, "S2C_priority_pair_pathways.png"), p_bubble, width = 10, height = 8, dpi = 300)
}

focus_source <- "Proliferative tumor"
focus_target <- "Immunosuppressive / M2-like TAM"

focus_lr <- lr_all %>%
  filter(source == focus_source, target == focus_target) %>%
  arrange(desc(prob)) %>%
  slice_head(n = 5)

focus_pathway <- pathway_all %>%
  filter(source == focus_source, target == focus_target) %>%
  arrange(desc(prob)) %>%
  slice_head(n = 5)

write.csv(focus_lr, file.path(round1_dir, "Fig2E_top5_LR_Prolif_to_M2TAM.csv"), row.names = FALSE)
write.csv(focus_pathway, file.path(round1_dir, "Fig2F_top5_pathways_Prolif_to_M2TAM.csv"), row.names = FALSE)

if (nrow(focus_lr) > 0) {
  focus_lr$interaction_name <- factor(
    focus_lr$interaction_name,
    levels = rev(focus_lr$interaction_name)
  )
  p_lr <- ggplot(
    focus_lr,
    aes(x = prob, y = interaction_name, size = prob, color = prob)
  ) +
    geom_point(alpha = 0.95) +
    scale_color_viridis_c(option = "D", name = "Probability") +
    scale_size_continuous(range = c(4, 10), name = "Probability") +
    labs(x = "Communication probability", y = "Ligand-receptor pair") +
    theme_classic(base_size = 14)
  ggsave(file.path(round1_dir, "Fig2E_top5_LR_Prolif_to_M2TAM.pdf"), p_lr, width = 8, height = 5)
}

if (nrow(focus_pathway) > 0) {
  focus_pathway$pathway_name <- factor(
    focus_pathway$pathway_name,
    levels = rev(focus_pathway$pathway_name)
  )
  p_pathway <- ggplot(
    focus_pathway,
    aes(x = prob, y = pathway_name, fill = prob)
  ) +
    geom_col(width = 0.75) +
    scale_fill_viridis_c(option = "D", name = "Strength") +
    labs(x = "Communication strength", y = "Pathway") +
    theme_classic(base_size = 14)
  ggsave(file.path(round1_dir, "Fig2F_top5_pathways_Prolif_to_M2TAM.pdf"), p_pathway, width = 8, height = 5)
}

# Round 2: remove T/NK cells and dendritic cells before recomputing CellChat.
immune_exclude <- c("T / NK cells", "Dendritic cells")
obj2 <- subset(obj, subset = !celltype_refined %in% immune_exclude)
obj2$celltype_refined <- droplevels(factor(obj2$celltype_refined))

cellchat2 <- run_cellchat(obj2)
saveRDS(cellchat2, file.path(round2_dir, "cellchat_no_TNK_DC.rds"))

write.csv(subsetCommunication(cellchat2), file.path(round2_dir, "CellChat_no_TNK_DC_LR.csv"), row.names = FALSE)
write.csv(subsetCommunication(cellchat2, slot.name = "netP"), file.path(round2_dir, "CellChat_no_TNK_DC_pathway.csv"), row.names = FALSE)
write.csv(cellchat2@net$count, file.path(round2_dir, "CellChat_no_TNK_DC_count_matrix.csv"))
write.csv(cellchat2@net$weight, file.path(round2_dir, "CellChat_no_TNK_DC_weight_matrix.csv"))

tumor_subtypes <- c(
  "Proliferative tumor",
  "Mesenchymal/Invasive tumor",
  "ECM-rich/claudin-high tumor",
  "Stress-hypoxia epithelial-like tumor",
  "Neural/Stem-like tumor"
)
myeloid_subtypes <- c(
  "Resident/homeostatic macrophage",
  "Immunosuppressive / M2-like TAM",
  "Activated TAM",
  "Metallothionein / APC-like macrophage",
  "Inflammatory myeloid",
  "IFN-responsive myeloid",
  "Neutrophil-like inflammatory myeloid",
  "Cycling myeloid"
)

groups <- colnames(cellchat2@net$weight)
tumor_present <- intersect(tumor_subtypes, groups)
myeloid_present <- intersect(myeloid_subtypes, groups)
weight <- cellchat2@net$weight

plot_strength_heatmap <- function(mat, x_name, y_name, title, filename) {
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("Sender", "Receiver", "Strength")
  p <- ggplot(df, aes(x = Receiver, y = Sender, fill = Strength)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(Strength, 3)), size = 3.5) +
    scale_fill_viridis_c(option = "D") +
    labs(title = title, x = x_name, y = y_name, fill = "Strength") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  ggsave(file.path(round2_dir, filename), p, width = 10, height = 8)
}

plot_strength_heatmap(
  weight[tumor_present, myeloid_present, drop = FALSE],
  "Myeloid subtype (receiver)",
  "Tumor subtype (sender)",
  "Tumor subtypes -> Myeloid subtypes",
  "Fig2D_tumor_to_myeloid_strength.pdf"
)

plot_strength_heatmap(
  weight[myeloid_present, tumor_present, drop = FALSE],
  "Tumor subtype (receiver)",
  "Myeloid subtype (sender)",
  "Myeloid subtypes -> Tumor subtypes",
  "S2D_myeloid_to_tumor_strength.pdf"
)
