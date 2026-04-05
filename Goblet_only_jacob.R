# =============================================================================
# GOBLET CELL SUB-ANNOTATION SCRIPT (Goblet-Only RDS — Jacob)
# =============================================================================
# This script loads the pre-processed goblet-only RDS (4.7 GB) which already
# contains UMAP coordinates and clusters. It skips the subsetting/re-clustering
# from the full 15.7 GB file and goes straight to annotation + plots.
# =============================================================================

# --- SETUP ---
library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(writexl)
library(tidyr)
library(ggpubr)

root_path <- "/Users/jacobng/Research/cai/wu_directory/wu_project"
output_dir <- file.path(root_path, "seurat_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(root_path)

# =============================================================================
# --- LOAD DATA (Goblet-only RDS) ---
# =============================================================================
message("Loading goblet-only RDS file...")
data_goblet <- readRDS(file.path(output_dir, "data_wu_project_goblet_subpops_final.rds"))
message(paste("Loaded object with", ncol(data_goblet), "cells and", nrow(data_goblet), "genes"))

# --- Discover what's in the object ---
message("\n--- Available reductions ---")
print(names(data_goblet@reductions))
message("\n--- Available metadata columns ---")
print(colnames(data_goblet@meta.data))
message("\n--- Existing cluster/subtype counts ---")
if ("seurat_clusters" %in% colnames(data_goblet@meta.data)) {
  message("seurat_clusters:")
  print(table(data_goblet$seurat_clusters))
}
if ("goblet_subtypes" %in% colnames(data_goblet@meta.data)) {
  message("goblet_subtypes (from colleague):")
  print(table(data_goblet$goblet_subtypes))
}

# --- Determine which UMAP reduction to use ---
# The colleague's object may have "umap", "umap_harmony", or both
umap_reduction <- if ("umap" %in% names(data_goblet@reductions)) {
  "umap"
} else if ("umap_harmony" %in% names(data_goblet@reductions)) {
  "umap_harmony"
} else {
  stop("No UMAP reduction found in the RDS! Available: ", paste(names(data_goblet@reductions), collapse = ", "))
}
message(paste("Using UMAP reduction:", umap_reduction))

# --- Determine which cluster column to use ---
cluster_col <- if ("seurat_clusters" %in% colnames(data_goblet@meta.data)) {
  "seurat_clusters"
} else if ("goblet_subtypes" %in% colnames(data_goblet@meta.data)) {
  "goblet_subtypes"
} else {
  stop("No cluster column found! Available columns: ", paste(colnames(data_goblet@meta.data), collapse = ", "))
}
message(paste("Using cluster column:", cluster_col))

# =============================================================================
# --- DEFINE FUNCTIONS ---
# =============================================================================

# 1. Weighted Annotation Function (same as Processing(Mine).R)
get_weighted_annotation <- function(seurat_obj, marker_genes, cluster_key, standardize_expression = TRUE) {
  all_obj_genes <- rownames(seurat_obj)
  marker_genes <- lapply(marker_genes, function(gs) intersect(gs, all_obj_genes))
  
  # Report which genes were found/missing per subtype
  for (ct in names(marker_genes)) {
    found <- marker_genes[[ct]]
    original <- marker_genes[[ct]]  # after intersect
    message(paste0("  ", ct, ": ", length(found), " genes found"))
  }
  
  marker_genes <- marker_genes[sapply(marker_genes, length) > 0]
  all_marker_genes <- sort(unique(unlist(marker_genes)))
  cell_types <- names(marker_genes)
  scores_binary <- matrix(0, nrow = length(all_marker_genes), ncol = length(cell_types), dimnames = list(all_marker_genes, cell_types))
  for (ct in cell_types) { scores_binary[marker_genes[[ct]], ct] <- 1.0 }
  W_specificity <- scores_binary / rowSums(scores_binary)
  W_final <- sweep(W_specificity, 2, colSums(W_specificity), "/")
  W_final[is.na(W_final)] <- 0
  X_subset <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[all_marker_genes, ]
  if (standardize_expression) { X_subset <- t(scale(t(as.matrix(X_subset)))); X_subset[is.na(X_subset)] <- 0 }
  clusters <- seurat_obj[[cluster_key, drop = TRUE]]
  annotation_vector <- character(length = ncol(seurat_obj)); names(annotation_vector) <- colnames(seurat_obj)
  for (cluster_id in unique(clusters)) {
    cell_idx <- which(clusters == cluster_id)
    mean_scores <- rowMeans(t(W_final) %*% X_subset[, cell_idx, drop = FALSE])
    annotation_vector[cell_idx] <- names(mean_scores)[which.max(mean_scores)]
  }
  return(list(annotation_vector = annotation_vector))
}

# 2. Gene Comparison Plot Function (same as Processing(Mine).R)
generate_gene_comparison_plots <- function(seurat_obj, score_col, group_by, x_axis, comparisons, plot_type = "violin", output_prefix, plot_title, y_label, fig_width = 12, fig_height = 8) {
  df_plot <- FetchData(seurat_obj, vars = c(score_col, group_by, x_axis))
  colnames(df_plot)[1] <- "Expression"
  p <- ggplot(df_plot, aes(x = !!sym(x_axis), y = Expression, fill = !!sym(x_axis))) +
    geom_violin(trim = TRUE, alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white") +
    stat_compare_means(comparisons = comparisons, label = "p.signif") +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_y") +
    labs(title = plot_title, y = y_label) + theme_classic()
  ggsave(file.path(output_dir, paste0(output_prefix, score_col, ".png")), p, width = fig_width, height = fig_height)
  return(p)
}

# =============================================================================
# --- PLOT 1: UMAP OF CLUSTERS (before annotation) ---
# =============================================================================
message("Generating UMAP of clusters...")
p_clusters <- DimPlot(data_goblet, reduction = umap_reduction, group.by = cluster_col, label = TRUE, repel = TRUE) +
  ggtitle("Goblet Cell Clusters") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
ggsave(file.path(output_dir, "UMAP_goblet_clusters.png"), p_clusters, width = 10, height = 7, dpi = 300, bg = "white")
message("  -> Saved: UMAP_goblet_clusters.png")

# =============================================================================
# --- PLOT 2: DOTPLOT OF ALL CANDIDATE MARKERS BY CLUSTER ---
# =============================================================================
message("Generating dotplot of all candidate markers...")

# All markers from the 8-subtype annotation guide (docx)
all_goblet_markers <- c(
  # Canonical Crypt GCs
  "Clca1", "Fcgbp", "Atoh1", "Spdef",
  # Noncanonical GCs
  "Dmbt1", "Gsdmc4", "Aqp8", "Muc17",
  # Intercrypt GCs (icGCs)
  "Slfn4", "Stxbp1", "Mxd1", "Ido1",
  # Proliferative GCs
  "Mki67", "Top2a", "Pcna", "Stmn1",
  # Sentinel GCs (senGCs)
  "Nlrp6", "Il18", "Wfdc2", "Areg",
  # RELMb+ GCs
  "Retnlb",
  # TFF3+ Repair GCs
  "Tff3", "Tff1",
  # MUC5AC+ Metaplastic GCs
  "Muc5ac", "Muc6",
  # Pan-GC reference
  "Muc2",
  # OLD "Defense" markers (Paneth cell) — included to SHOW they have no signal
  "Lyz1", "Defa5", "Defa6"
)

# Filter to genes actually present in dataset
genes_present <- intersect(all_goblet_markers, rownames(data_goblet))
genes_missing <- setdiff(all_goblet_markers, rownames(data_goblet))
message(paste("  Markers found:", length(genes_present), "/ Missing:", length(genes_missing)))
if (length(genes_missing) > 0) {
  message(paste("  Missing genes:", paste(genes_missing, collapse = ", ")))
}

dot_plot <- DotPlot(object = data_goblet, dot.min = 0.05, cols = "RdBu",
                    features = genes_present, scale = TRUE,
                    group.by = cluster_col) + coord_flip() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    strip.text = element_blank()
  ) +
  ggtitle("Goblet Sub-type Markers by Cluster")

ggsave(file.path(output_dir, "Dotplot_goblet_subtype_markers_by_cluster.png"),
       plot = dot_plot, width = 10, height = 14, dpi = 300, bg = "white")
message("  -> Saved: Dotplot_goblet_subtype_markers_by_cluster.png")

# =============================================================================
# --- GOBLET ANNOTATION WITH CORRECTED 8-SUBTYPE MARKERS ---
# =============================================================================
message("Running weighted annotation with corrected 8-subtype markers...")

goblet_sub_markers <- list(
  "Canonical"       = c("Clca1", "Fcgbp", "Atoh1", "Spdef"),
  "Noncanonical"    = c("Dmbt1", "Gsdmc4", "Aqp8", "Muc17"),
  "icGC"            = c("Slfn4", "Stxbp1", "Mxd1", "Ido1"),
  "Proliferative"   = c("Mki67", "Top2a", "Pcna", "Stmn1"),
  "Sentinel"        = c("Nlrp6", "Il18", "Wfdc2", "Areg"),
  "RELMb+"          = c("Retnlb", "Clca1"),
  "TFF3+ Repair"    = c("Tff3", "Tff1"),
  "MUC5AC+ Meta"    = c("Muc5ac", "Muc6")
)

results <- get_weighted_annotation(data_goblet, goblet_sub_markers, cluster_col)
data_goblet$sub_cell_types <- results$annotation_vector

# Print summary
message("\n--- Annotation Summary ---")
print(table(data_goblet$sub_cell_types))

# =============================================================================
# --- PLOT 3: UMAP WITH SUBTYPE ANNOTATIONS ---
# =============================================================================
message("Generating annotated UMAP...")
p_annotated <- DimPlot(data_goblet, reduction = umap_reduction, group.by = "sub_cell_types", label = TRUE, repel = TRUE) +
  ggtitle("Goblet Sub-populations (Corrected Annotation)") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
ggsave(file.path(output_dir, "UMAP_goblet_sub_types.png"), p_annotated, width = 10, height = 7, dpi = 300, bg = "white")
message("  -> Saved: UMAP_goblet_sub_types.png")

# =============================================================================
# --- PLOT 3b: UMAP SPLIT BY GENOTYPE x DIET (2x2 grid) ---
# =============================================================================
message("Generating 2x2 UMAP split by Genotype_Diet...")
p_split <- DimPlot(
  object = data_goblet,
  reduction = umap_reduction,
  group.by = "sub_cell_types",
  split.by = "Genotype_Diet",
  label = FALSE,
  repel = TRUE
) +
  facet_wrap(~Genotype_Diet, ncol = 2) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(file.path(output_dir, "UMAP_goblet_sub_types_2x2.png"),
       plot = p_split, width = 12, height = 10, dpi = 300, bg = "white")
message("  -> Saved: UMAP_goblet_sub_types_2x2.png")

# =============================================================================
# --- PLOT 4: DOTPLOT OF MARKERS BY ANNOTATED SUBTYPE ---
# =============================================================================
message("Generating dotplot by annotated subtype...")
dot_plot_annotated <- DotPlot(object = data_goblet, dot.min = 0.05, cols = "RdBu",
                               features = genes_present[!genes_present %in% c("Lyz1", "Defa5", "Defa6")],
                               scale = TRUE,
                               group.by = "sub_cell_types") + coord_flip() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    strip.text = element_blank()
  ) +
  ggtitle("Goblet Sub-type Markers by Annotation")

ggsave(file.path(output_dir, "Dotplot_goblet_markers_by_subtype.png"),
       plot = dot_plot_annotated, width = 10, height = 12, dpi = 300, bg = "white")
message("  -> Saved: Dotplot_goblet_markers_by_subtype.png")

# =============================================================================
# --- PLOT 5: PPARG COMPARISON ---
# =============================================================================
message("Generating Pparg comparison plots...")
p_pparg <- generate_gene_comparison_plots(
  data_goblet, "Pparg", "sub_cell_types", "Genotype_Diet",
  list(c("WT_inulin", "KO_inulin"), c("WT_cellulose", "KO_cellulose")),
  output_prefix = "Goblet_Sub_",
  plot_title = "Pparg Expression across Goblet Sub-populations",
  y_label = "Log-Normalized Expression",
  fig_width = 14, fig_height = 10
)
message("  -> Saved: Goblet_Sub_Pparg.png")

# =============================================================================
# --- SAVE ANNOTATED OBJECT ---
# =============================================================================
saveRDS(data_goblet, file.path(output_dir, "data_wu_project_goblet_sub_annotated.rds"))
message("\n=== Goblet analysis complete. All plots saved to seurat_output/ ===")
