# --- SETUP ---
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(writexl)
library(tidyr)
library(ggpubr)

root_path <- "/Users/aayush/RStudio"
output_dir <- file.path(root_path, "seurat_output")
setwd(root_path)

# --- LOAD DATA ---
# Using the correct filename based on your previous logs
message("Loading data...")
data <- readRDS(file.path(output_dir, "DecontX_Doubletfinders_Colonocytes.rds"))

# --- DEFINE MISSING FUNCTIONS ---

# 1. Weighted Annotation Function
get_weighted_annotation <- function(seurat_obj, marker_genes, cluster_key, standardize_expression = TRUE) {
  all_obj_genes <- rownames(seurat_obj)
  marker_genes <- lapply(marker_genes, function(gs) intersect(gs, all_obj_genes))
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

# 2. Gene Comparison Plot Function
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

# --- GOBLET SUB-CLUSTERING ---

process_and_extract_cell_type <- function(data_obj, cell_type_name) {
  # Change param name to 'data_obj' to avoid confusion with the 'data' function
  message(paste("Subsetting for:", cell_type_name))
  # Use broad_cell_types (weighted or original) as identified in previous step
  data_subset <- subset(data_obj, subset = broad_cell_types == cell_type_name)
  
  data_subset <- FindVariableFeatures(data_subset, nfeatures = 2000)
  data_subset <- ScaleData(data_subset)
  data_subset <- RunPCA(data_subset, npcs = 30, verbose = FALSE)
  data_subset <- RunHarmony(data_subset, group.by.vars = "SampleID", dims.use = 1:30, verbose = FALSE)
  data_subset <- FindNeighbors(data_subset, reduction = "harmony", dims = 1:30)
  data_subset <- FindClusters(data_subset, resolution = 0.8)
  data_subset <- RunUMAP(data_subset, reduction = "harmony", dims = 1:30)
  return(data_subset)
}

# Execute processing
data_goblet <- process_and_extract_cell_type(data, "Goblet cells")

# --- GOBLET ANNOTATION ---

goblet_sub_markers <- list(
  "Canonical" = c("Atoh1", "Muc2", "Fcgbp", "Clca1"),
  "Non-canonical" = c("Hes1", "Dmbt1", "Muc17", "Slc12a2"),
  "Cycling" = c("Mki67", "Top2a", "Pcna", "Stmn1"),
  "Defense" = c("Lyz1", "Defa5", "Defa6", "Ccn3")
)

results <- get_weighted_annotation(data_goblet, goblet_sub_markers, "seurat_clusters")
data_goblet$sub_cell_types <- results$annotation_vector

# --- FINAL PLOTS ---

# 1. UMAP
p1 <- DimPlot(data_goblet, reduction = "umap", group.by = "sub_cell_types", label = TRUE) + ggtitle("Goblet Sub-populations")
ggsave(file.path(output_dir, "UMAP_goblet_sub_types.png"), p1, width = 10, height = 7)

# 2. Pparg Comparison
p2 <- generate_gene_comparison_plots(
  data_goblet, "Pparg", "sub_cell_types", "Genotype_Diet",
  list(c("WT_inulin", "KO_inulin"), c("WT_cellulose", "KO_cellulose")),
  output_prefix = "Goblet_Sub_"
)

# 3. Save
saveRDS(data_goblet, file.path(output_dir, "data_wu_project_goblet_sub_annotated.rds"))
message("Goblet analysis complete.")