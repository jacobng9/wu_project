# =============================================================================
# MASTER SCRIPT FOR SINGLE-CELL ANALYSIS (WU IEC PROJECT) -Annotation
# =============================================================================

# --- 1. SETUP: LOAD LIBRARIES AND DEFINE PATHS ---
library(Seurat)
library(harmony)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(patchwork)
library(writexl)
library(celda)
library(DoubletFinder)
library(Matrix)

set.seed(123)
root_path <- "Z:/selim_working_dir/2025_wu_iec_project/r_process"
#root_path <- "/mnt/SCDC/Optimus/selim_working_dir/2025_wu_iec_project/r_process" 
setwd(root_path)
output_dir <- file.path(root_path, "seurat_output")

data<- readRDS(file.path(output_dir, "data_wu_project_decontX_doubletfinders.rds"))

# To see the actual clusters side-by-side
p1 <- DimPlot(data, reduction = "umap_harmony", group.by = "clusters_harmony", label = TRUE)
p1
ggsave(file.path(output_dir, "UMAP_Cluster_for_annotation.png"), plot = p1, width = 8, height = 7, dpi = 300)
# To see the actual clusters side-by-side
p1 <- DimPlot(data, reduction = "umap_harmony", group.by = "SampleID", label = FALSE)
p1
ggsave(file.path(output_dir, "UMAP_SampleID_comps.png"), plot = p1, width = 8, height = 7, dpi = 300)


cell_type_markers <- c(
  "Ighg1", "Mzb1", # Plasma B cells
  "Ms4a1", "Cd19", # B cells
  "Cd3e", "Cd4", "Cd8a", "Trbc2", # T cells
  "Clec9a", "Xcr1", # cDCs
  "Cd163", "Csf1r", "Cd68", "Mrc1", "Aif1", # Macrophages
  "Col1a1", "Col3a1", "Col6a2", # Fibroblasts
  "Vip", "Tubb3", "Vat1l", # ENs   - Enteric neurons
  "Sox10", "Plp1", # EGCs    - Enteric glial cells
  "Lyve1", "Flt4", "Pecam1", # LECs      -   Lymphatic endothelial cells
  "Plvap", "Flt1", # VECs    -   Vascular endothelial cells
  "Epcam", "Muc2", "Krt19", "Krt20", "Vil1", # Colonocytes
  "Myh11", "Tagln", "Hhip", "Des", # SMCs    -   Smooth muscle cells
  "Adipoq", "Plin1", "Plin4" # Adipocytes
)

# Create the DotPlot
dot_plot <- DotPlot(object = data, dot.min = 0.05, cols = 'RdBu',
                    features = cell_type_markers, scale = T,
                    group.by = "clusters_harmony") +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 13),
    axis.text.y = element_text(size = 15),
    strip.text = element_blank()
  ) 
# Display the plot
dot_plot
ggsave(file.path(output_dir, "Marker_dotplot_harmony_clusters.png"), 
       plot = dot_plot, width = 11, height = 12, dpi = 300, bg = "white")


library(dplyr)
library(tidyr)
library(Matrix)

get_weighted_annotation <- function(
    seurat_obj, 
    marker_genes, 
    cluster_key, 
    standardize_expression = TRUE
) {
  # --- Step 0: Sync markers with dataset ---
  all_obj_genes <- rownames(seurat_obj)
  marker_genes <- lapply(marker_genes, function(gs) intersect(gs, all_obj_genes))
  marker_genes <- marker_genes[sapply(marker_genes, length) > 0]
  
  all_marker_genes <- sort(unique(unlist(marker_genes)))
  cell_types <- names(marker_genes)
  
  # --- Step 1: Weighted Score Matrix (Identical to Python Logic) ---
  scores_binary <- matrix(0, nrow = length(all_marker_genes), ncol = length(cell_types),
                          dimnames = list(all_marker_genes, cell_types))
  for (ct in cell_types) { scores_binary[marker_genes[[ct]], ct] <- 1.0 }
  
  # Specificity weighting: 1 / gene_occurrence
  gene_occurrence <- rowSums(scores_binary)
  W_specificity <- scores_binary / gene_occurrence
  
  # Marker set normalization: Ensure column sums = 1
  column_sums <- colSums(W_specificity)
  W_final <- sweep(W_specificity, 2, column_sums, "/")
  W_final[is.na(W_final)] <- 0
  
  # --- Step 2: Standardize Expression (Z-score per gene across cells) ---
  # Pulling from 'data' layer (normalized counts)
  X_subset <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[all_marker_genes, ]
  
  if (standardize_expression) {
    # Scale works on columns, so we transpose twice to Z-score genes
    X_subset <- t(scale(t(as.matrix(X_subset)))) 
    X_subset[is.na(X_subset)] <- 0
  } else {
    X_subset <- as.matrix(X_subset)
  }
  
  # --- Step 3: Scoring per Cluster ---
  W <- t(W_final) # Result: (CellType x Gene)
  clusters <- seurat_obj[[cluster_key, drop = TRUE]]
  unique_clusters <- unique(clusters)
  
  # Initialize the final cell-level vector and a list for the ranking report
  annotation_vector <- character(length = ncol(seurat_obj))
  names(annotation_vector) <- colnames(seurat_obj)
  cluster_score_data <- list()
  
  for (cluster_id in unique_clusters) {
    cell_idx <- which(clusters == cluster_id)
    if(length(cell_idx) == 0) next
    
    # (CellType x Gene) @ (Gene x Cells_in_cluster)
    raw_scores <- W %*% X_subset[, cell_idx, drop = FALSE]
    mean_scores <- rowMeans(raw_scores)
    
    # Sort to find Top 5 (matches Python np.argsort()[::-1])
    top_indices <- order(mean_scores, decreasing = TRUE)
    top_ct <- names(mean_scores)[top_indices]
    top_sc <- mean_scores[top_indices]
    
    # Assign the #1 hit to ALL cells in this cluster (the annotation_vector)
    annotation_vector[cell_idx] <- top_ct[1]
    
    # Store ranking data
    n_top <- min(5, length(top_ct))
    cluster_score_data[[as.character(cluster_id)]] <- data.frame(
      cluster = cluster_id,
      Rank = 1:n_top,
      Cell_Type = top_ct[1:n_top],
      Score = top_sc[1:n_top]
    )
  }
  
  top5_scores_df <- do.call(rbind, cluster_score_data)
  
  # Return ONLY the vector and the report (no Seurat object duplication)
  return(list(
    annotation_vector = annotation_vector,
    top5_report = top5_scores_df
  ))
}

summarize_annotations <- function(seurat_obj, annotation_column, print_summary = TRUE) {
  if (!annotation_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", annotation_column, "not found in metadata."))
  }
  
  counts <- table(seurat_obj[[annotation_column]])
  percentages <- prop.table(counts) * 100
  
  summary_df <- data.frame(
    Count = as.numeric(counts),
    Percentage = as.numeric(percentages),
    row.names = names(counts)
  ) %>% arrange(desc(Percentage))
  
  if (print_summary) {
    print_df <- summary_df
    print_df$Percentage <- sprintf("%.2f%%", print_df$Percentage)
    message(paste("\n--- Annotation Summary for '", annotation_column, "' ---"))
    print(print_df)
  }
  
  return(summary_df)
}

# Define your marker list
my_markers <- list(
  "Plasma B cells" = c("Ighg1", "Mzb1"),
  "B cells"        = c("Ms4a1", "Cd19"),
  "T cells"        = c("Cd3e", "Cd4", "Cd8a", "Trbc2"),
  "cDCs"           = c("Clec9a", "Xcr1"),
  "Macrophages"    = c("Cd163", "Csf1r", "Cd68", "Mrc1", "Aif1"),
  "Fibroblasts"    = c("Col1a1", "Col3a1", "Col6a2"),
  "ENs"            = c("Vip", "Tubb3", "Vat1l"),
  "EGCs"           = c("Sox10", "Plp1"),
  "LECs"           = c("Lyve1", "Flt4", "Pecam1"),
  "VECs"           = c("Plvap", "Flt1"),
  "Colonocytes"    = c("Epcam", "Muc2", "Krt19", "Krt20", "Vil1"),
  "SMCs"           = c("Myh11", "Tagln", "Hhip", "Des"),
  "Adipocytes"     = c("Adipoq", "Plin1", "Plin4")
)

# 1. Run the function
# It will return a list with the per-cell vector and the ranking report
results <- get_weighted_annotation(
  seurat_obj = data, 
  marker_genes = my_markers, 
  cluster_key = "clusters_harmony",
  standardize_expression = TRUE
)

# 2. Assign the vector directly to metadata
# This vector is already ordered by cell, so it fits perfectly
data$broad_cell_types_weighted <- results$annotation_vector

# 3. Check the ranking report to see the scores for each cluster
print(head(results$top5_report))

res_sum_ann <- summarize_annotations(data, "broad_cell_types_weighted",print_summary = TRUE)
res_sum_ann
# 4. Final Visualization
DimPlot(data, reduction = "umap_harmony", group.by = "broad_cell_types_weighted", label = TRUE)



# 1. Run the function
# It will return a list with the per-cell vector and the ranking report
results <- get_weighted_annotation(
  seurat_obj = data, 
  marker_genes = my_markers, 
  cluster_key = "clusters_harmony",
  standardize_expression = FALSE
)

# 2. Assign the vector directly to metadata
# This vector is already ordered by cell, so it fits perfectly
data$broad_cell_types_weighted_not_scaled <- results$annotation_vector

# 3. Check the ranking report to see the scores for each cluster
print(head(results$top5_report))

res_sum_ann <- summarize_annotations(data, "broad_cell_types_weighted_not_scaled",print_summary = TRUE)
res_sum_ann
# 4. Final Visualization
DimPlot(data, reduction = "umap_harmony", group.by = "broad_cell_types_weighted_not_scaled", label = TRUE)




# Re-annotation
data$broad_cell_types = data$clusters_harmony
data$broad_cell_types = recode_factor(data$broad_cell_types,
                                      '0' = 'Colonocytes',
                                      '1' = 'Colonocytes',
                                      '2' = 'Colonocytes', 
                                      '3' = 'Colonocytes',
                                      '4' = 'Colonocytes',
                                      '5' = 'Colonocytes',
                                      '6' = 'Colonocytes', 
                                      '7' = 'Colonocytes',
                                      '8' = 'Plasma B cells',
                                      '9' = 'Colonocytes',
                                      '10' = 'Colonocytes',
                                      '11' = 'B cells',
                                      '12' = 'Colonocytes',
                                      '13' = 'Colonocytes',
                                      '14' = 'T cells',
                                      '15' = 'Macrophages',
                                      '16' = 'Fibroblasts',
                                      '17' = 'SMCs',
                                      '18' = 'Colonocytes',
                                      '19' = 'T cells',
                                      '20' = 'SMCs',
                                      '21' = 'Fibroblasts',
                                      '22' = 'VECs',
                                      '23' = 'cDCs',
                                      '24' = 'Colonocytes',
                                      '25' = 'SMCs',
                                      '26' = 'LECs',
                                      '27' = 'Plasma B cells',# or B cells
                                      '28' = 'Colonocytes',
                                      '29' = 'SMCs', 
                                      '30' = 'Plasma B cells',
                                      '31' = 'SMCs',
                                      '32' = 'SMCs',
                                      '33' = 'Colonocytes',
                                      '34' = 'Colonocytes',
                                      '35' = 'Colonocytes',
                                      '36' = 'Colonocytes',
                                      '37' = 'Colonocytes',
                                      '38' = 'EGCs',
                                      '39' = 'B cells',
                                      '40' = 'B cells',
                                      '41' = 'Fibroblasts',
                                      '42' = 'ENs',
                                      '43' = 'Colonocytes',
                                      '44' = 'Adipocytes',
                                      '45' = 'LECs',
                                      '46' = 'Colonocytes',
                                      '47' = 'Colonocytes',
                                      '48' = 'Colonocytes',
                                      '49' = 'Colonocytes',
                                      '50' = 'SMCs',
                                      '51' = 'LECs',
                                      '52' = 'T cells',
                                      '53' = 'Colonocytes'
                                      
)

DimPlot(object = data, reduction = "umap_harmony", group.by = c("broad_cell_types"), label = TRUE, repel = TRUE)

cell_type_markers <- c(
  "Ighg1", "Mzb1", # Plasma B cells
  "Ms4a1", "Cd19", # B cells
  "Cd3e", "Cd4", "Cd8a", "Trbc2", # T cells
  "Clec9a", "Xcr1", # cDCs
  "Cd163", "Csf1r", "Cd68", "Mrc1", "Aif1", # Macrophages
  "Col1a1", "Col3a1", "Col6a2", # Fibroblasts
  "Vip", "Tubb3", "Vat1l", # ENs   - Enteric neurons
  "Sox10", "Plp1", # EGCs    - Enteric glial cells
  "Lyve1", "Flt4", "Pecam1", # LECs      -   Lymphatic endothelial cells
  "Plvap", "Flt1", # VECs    -   Vascular endothelial cells
  "Epcam", "Muc2", "Krt19", "Krt20", "Vil1", # Colonocytes
  "Myh11", "Tagln", "Hhip", "Des", # SMCs    -   Smooth muscle cells
  "Adipoq", "Plin1", "Plin4" # Adipocytes
)

cell_levels <- c("Plasma B cells", "B cells", "T cells", "cDCs", "Macrophages", 
                 "Fibroblasts", "ENs", "EGCs", "LECs", "VECs", "Colonocytes", 
                 "SMCs", "SMCs2", "Adipocytes")

data$broad_cell_types <- factor(data$broad_cell_types, levels = cell_levels)
data$broad_cell_types_weighted <- factor(data$broad_cell_types_weighted, levels = cell_levels)
data$broad_cell_types_weighted_not_scaled <- factor(data$broad_cell_types_weighted_not_scaled, levels = cell_levels)

# UMAP of broad cell types
p1 <- DimPlot(object = data, reduction = "umap_harmony", 
              group.by = c("broad_cell_types",
                           "broad_cell_types_weighted",
                           "broad_cell_types_weighted_not_scaled"), 
              label = TRUE, repel = TRUE) +
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16, face = "bold")
  ) + guides(color = guide_legend(override.aes = list(size = 4)))
p1 

# UMAP of broad cell types assigning winner
p1 <- DimPlot(object = data, reduction = "umap_harmony", 
              group.by = c("broad_cell_types"), 
              label = TRUE, repel = TRUE) +
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16, face = "bold")
  ) + guides(color = guide_legend(override.aes = list(size = 4)))
p1 
ggsave(file.path(output_dir, "UMAP_broad_cell_types.png"), 
       plot = p1, width = 10, height = 7, dpi = 300, bg = "white")



# Create the DotPlot of broad cell type gene markers
dot_plot <- DotPlot(object = data, dot.min = 0.05, cols = "RdBu", #c("white","red"),
                    features = cell_type_markers, scale = T,
                    group.by = "broad_cell_types") + coord_flip() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 13),
    axis.text.y = element_text(size = 15),
    strip.text = element_blank()
  ) 
# Display the plot
dot_plot
ggsave(file.path(output_dir, "Dotplot_markers_broad_cell_types.png"), 
       plot = dot_plot, width = 7, height = 10, dpi = 300, bg = "white")



# Create the vector and assign it directly having Genotype_Diet effect 
library(dplyr)
data$Genotype_Diet <- ifelse(grepl("vWT", data$Batch) & data$Diet == "cellulose", "WT_cellulose",
                             ifelse(grepl("vKO", data$Batch) & data$Diet == "cellulose", "KO_cellulose",
                                    ifelse(grepl("vWT", data$Batch) & data$Diet == "inulin", "WT_inulin",
                                           ifelse(grepl("vKO", data$Batch) & data$Diet == "inulin", "KO_inulin", "Unknown"))))

data$Genotype_Diet <- factor(data$Genotype_Diet, 
                             levels = c("WT_cellulose", "KO_cellulose", "WT_inulin", "KO_inulin"))

p1 <- DimPlot(
  object = data, 
  reduction = "umap_harmony", 
  group.by = "broad_cell_types", 
  split.by = "Genotype_Diet",  # This "injects" the variable into the plot layers
  label = FALSE, 
  repel = TRUE
) + 
  # 2. Now facet_wrap will recognize the variable
  facet_wrap(~Genotype_Diet, ncol = 2) + 
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16, face = "bold"),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(size = 14, face = "bold") 
  ) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

# Display
p1
# Update dimensions: For a 2x2 grid, a more square aspect ratio works best
ggsave(file.path(output_dir, "UMAP_broad_cell_types_2x2.png"), 
       plot = p1, width = 12, height = 10, dpi = 300, bg = "white")



library(dplyr)
library(tidyr)
library(ggpubr)            # Add this to the top of your script
library(ggplot2)
library(writexl) # Required for Excel export
# --- 1. THEME DEFINITION ---
custom_plot_theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom"
  )
plot_cell_proportions <- function(seurat_obj, cluster_col, output_prefix, output_dir) {
  
  # 1. Prepare Data
  metadata <- seurat_obj@meta.data
  make_label <- function(p) { ifelse(p > 2.5, paste0(round(p, 1), "%"), "") }
  
  # --- A. Stats by SampleID ---
  df_sample <- metadata %>%
    group_by(SampleID, !!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(SampleID) %>%
    mutate(percentage = (n / sum(n)) * 100)
  
  # --- B. Stats by Genotype_Diet ---
  df_group <- metadata %>%
    group_by(Genotype_Diet, !!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Genotype_Diet) %>%
    mutate(percentage = (n / sum(n)) * 100)
  
  # --- C. Global Stats ---
  df_global <- metadata %>%
    group_by(!!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(percentage = (n / sum(n)) * 100)
  
  # 2. Save Data to Excel (Multi-sheet)
  excel_path <- file.path(output_dir, paste0(output_prefix, "_Stats.xlsx"))
  
  # Create a named list where names = Sheet Names
  export_list <- list(
    "By_Sample" = df_sample,
    "By_Group"  = df_group,
    "Global"    = df_global
  )
  
  write_xlsx(export_list, path = excel_path)
  message("Excel saved: ", excel_path)
  
  # 3. Generate Plots (Keeping your previous logic)
  p1 <- ggplot(df_sample, aes(x = SampleID, y = percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
    geom_text(aes(label = make_label(percentage)), position = position_stack(vjust = 0.5), size = 4) +
    labs(title = paste("Proportions by Sample (", cluster_col, ")"), y = "Percentage (%)") +
    custom_plot_theme
  
  p2 <- ggplot(df_group, aes(x = Genotype_Diet, y = percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
    geom_text(aes(label = make_label(percentage)), position = position_stack(vjust = 0.5), size = 4) +
    labs(title = paste("Proportions by Group (", cluster_col, ")"), y = "Percentage (%)") +
    custom_plot_theme
  
  p3 <- ggplot(df_global, aes(x = "Global", y = percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", width = 0.6) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5), size = 5) +
    labs(title = paste("Global Distribution (", cluster_col, ")"), y = "Percentage (%)") +
    custom_plot_theme
  
  # 4. Save Plots
  ggsave(file.path(output_dir, paste0(output_prefix, "_cell_pct_by_Sample.png")), p1, width = 12, height = 8, bg = "white")
  ggsave(file.path(output_dir, paste0(output_prefix, "_cell_pct_by_Group.png")),  p2, width = 10, height = 8, bg = "white")
  ggsave(file.path(output_dir, paste0(output_prefix, "_cell_pct_Global.png")),    p3, width = 8, height = 8, bg = "white")
  
  return(list(plots = list(p1, p2, p3), stats = export_list))
}


# Define your labels to test
labels_to_test <- c("broad_cell_types")

# Loop through and generate all plots
for (label in labels_to_test) {
  message("Processing proportions for: ", label)
  plot_cell_proportions(
    seurat_obj = data, 
    cluster_col = label, 
    output_prefix = paste0("Colonocytes_", label), 
    output_dir = output_dir
  )
}


generate_gene_comparison_plots <- function(seurat_obj,
                                           score_col,
                                           group_by,
                                           x_axis,
                                           comparisons,
                                           plot_type = "violin",
                                           output_prefix,
                                           plot_title,
                                           y_label,
                                           fig_width = 16,  # New variable
                                           fig_height = 16  # New variable
) {
  
  message(paste("\n--- Processing:", score_col, "---"))
  
  # 1. Pull Data
  df_plot <- FetchData(seurat_obj, vars = c(score_col, group_by, x_axis))
  colnames(df_plot)[1] <- "Expression"
  df_plot <- df_plot %>% drop_na()
  
  # 2. Theme setup
  cust_theme <- theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 10)),
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1.2),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black", face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      panel.spacing = unit(1.5, "lines") 
    )
  
  # 3. Plotting Logic
  p <- ggplot(df_plot, aes(x = !!sym(x_axis), y = Expression, fill = !!sym(x_axis)))
  
  if (plot_type == "barplot") {
    p <- p +
      stat_summary(fun = mean, geom = "bar", color = "black", alpha = 0.8, linewidth = 0.8) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 0.8)
  } else {
    p <- p +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.7, linewidth = 0.8) +
      geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5, linewidth = 0.6)
  }
  
  # Add Stats and Formatting
  p <- p +
    stat_compare_means(
      comparisons = comparisons, 
      label = "p.signif", 
      method = "wilcox.test",
      method.args = list(exact = FALSE),
      # Standardize 'ns' labels
      symnum.args = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
        symbols = c("****", "***", "**", "*", "ns")
      ),
      step.increase = 0.1,  
      label.y.npc = "top", 
      size = 7,             
      bracket.size = 0.8    
    )  +
    facet_wrap(as.formula(paste("~", group_by)), scales = "free_y", ncol = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
    coord_cartesian(clip = "off") +
    labs(title = plot_title, y = y_label, fill = "Condition:") +
    scale_fill_brewer(palette = "Set1") +
    cust_theme +
    guides(fill = guide_legend(
      override.aes = list(alpha = 1), 
      keywidth = unit(2.5, "lines"), 
      keyheight = unit(1.5, "lines")
    ))
  
  # 4. Save Plot using the new variables
  if(!dir.exists("plots")) dir.create("plots")
  plot_path <- file.path("plots", paste0(output_prefix, score_col, ".png"))
  
  ggsave(plot_path, p, width = fig_width, height = fig_height, dpi = 300) 
  
  return(p)
}

# Example for a specific marker gene
p1<- generate_gene_comparison_plots(
  seurat_obj = data,
  score_col = "Pparg",
  group_by = "broad_cell_types",
  x_axis = "Genotype_Diet",
  comparisons = list(c("WT_inulin", "KO_inulin"), 
                     c("WT_cellulose", "KO_cellulose"), 
                     c("WT_cellulose", "WT_inulin")),
  plot_type = "violin",
  output_prefix = "pparg_mean_exp_bar_plots_broad_cell_types_",
  plot_title = "Pparg Expression across Groups",
  y_label = "Log-Normalized Expression"
)

p1




saveRDS(data, file.path(output_dir, "data_wu_project_decX_dblt_ann.rds"))
#saveRDS(data, file.path(output_dir, "data_wu_project_decX_dblt_ann.rds"), compress = TRUE)







######################################################
############# Colonocytes subpopulation annotations
#####################################################
library(Seurat)
library(harmony)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(patchwork)
library(writexl)
set.seed(123)
root_path <- "Z:/selim_working_dir/2025_wu_iec_project/r_process"
setwd(root_path)
output_dir <- file.path(root_path, "seurat_output")

data<- readRDS(file.path(output_dir, "data_wu_project_decX_dblt_ann.rds"))

process_and_extract_cell_type <- function(data, cell_type_name, assay_name = "RNA", num_hvg = 5000, 
                                          dims_pca = 50, resolution = 1.0, min_dist= 0.2, kneigh=15) {
  DefaultAssay(data) <- assay_name
  
  message(paste("Subsetting for:", cell_type_name))
  # 1. Create the subset
  data_subset <- subset(data, subset = CellType == cell_type_name)
  
  # 2. Clean up old dimensionality reductions/graphs to avoid the warning/bloat
  data_subset@reductions <- list() 
  data_subset@graphs <- list()
  
  # 3. Re-run normalization/scaling on the subset
  # It is best practice to re-find HVGs for a specific sub-population
  data_subset <- FindVariableFeatures(data_subset, selection.method = "vst", nfeatures = num_hvg)
  data_subset <- ScaleData(data_subset) # This overwrites scale.data for the new HVGs
  data_subset <- RunPCA(data_subset, npcs = dims_pca, reduction.name = "pca")
  gc()
  # --- VERSION 1: NON-HARMONY ---
  # FIX: Use 'data_subset' here, not 'data'
  data_subset <- FindNeighbors(data_subset, dims = 1:dims_pca, reduction = "pca", k.param = kneigh, graph.name = "pca_nn")
  data_subset <- FindClusters(data_subset, resolution = resolution, graph.name = "pca_nn", cluster.name = "clusters_none")
  data_subset <- RunUMAP(data_subset, dims = 1:dims_pca, reduction = "pca", n.neighbors = kneigh, 
                         reduction.name = "umap_none", min.dist = min_dist, n.epochs = 500)
  gc()
  # --- VERSION 2: HARMONY ---
  message("Processing Harmony Pipeline...")
  library(harmony)
  data_subset <- RunHarmony(data_subset, group.by.vars = "SampleID", reduction = "pca", reduction.save = "harmony")
  gc()
  data_subset <- FindNeighbors(data_subset, dims = 1:dims_pca, reduction = "harmony", k.param = kneigh, graph.name = "harmony_nn")
  data_subset <- FindClusters(data_subset, resolution = resolution, graph.name = "harmony_nn", cluster.name = "clusters_harmony")
  data_subset <- RunUMAP(data_subset, dims = 1:dims_pca, reduction = "harmony", n.neighbors = kneigh, 
                         reduction.name = "umap_harmony", min.dist = min_dist, n.epochs = 500)
  gc()
  return(data_subset)
}

data$CellType <- data$broad_cell_types
data_colo <- process_and_extract_cell_type(data, "Colonocytes", num_hvg = 2000, 
                                           resolution = 3.0, min_dist = 0.3, kneigh = 15 )
# --- VISUAL COMPARISON ---
# Side-by-side comparison of Harmony vs. No Harmony
p1 <- DimPlot(data_colo, reduction = "umap_none", group.by = "SampleID") + ggtitle("Standard PCA (No Harmony)")
p2 <- DimPlot(data_colo, reduction = "umap_harmony", group.by = "SampleID") + ggtitle("Harmony Corrected")

# Combine plots and save
p_comp <- p1 + p2 & theme(legend.position = "bottom")
p_comp
ggsave(file.path(output_dir, "UMAP_Harmony_Comparison_colonocytes.png"), plot = p_comp, width = 16, height = 8, dpi = 300)


# To see the actual clusters side-by-side
p3 <- DimPlot(data_colo, reduction = "umap_none", group.by = "clusters_none", label = TRUE) + ggtitle("Clusters: Standard")
p4 <- DimPlot(data_colo, reduction = "umap_harmony", group.by = "clusters_harmony", label = TRUE) + ggtitle("Clusters: Harmony")

p_clusters <- p3 + p4
p_clusters
ggsave(file.path(output_dir, "UMAP_Cluster_Comparison_colonocytes.png"), plot = p_clusters, width = 16, height = 8, dpi = 300)




# Define your marker list
sub_markers <- list(
  "Abs. colonocytes"  = c("Alpi", "Ces2c", "Slc26a2", "Ceacam1", "Aqp8"),
  "Goblet cells"     = c("Muc2", "Spink4", "Agr2", "Tff3"),
  "EECs"             = c("Chga", "Chgb", "Tph1", "Cck", "Scgn", "Scg5"),
  "Tuft cells"       = c("Dclk1", "Trpm5", "Avil", "Sh2d6", "Plcg2"),
  "TA cells"         = c("Mki67", "Top2a", "Birc5", "Pcna", "Stmn1"),
  "Stem cells"       = c("Lgr5", "Lrig1", "Ascl2", "Slc12a2", "Smoc2", "Kcnq1", "Gpx2", "Ephb2", "Bmpr1a", "Hopx", "Sox9")
)

results <- get_weighted_annotation(
  seurat_obj = data_colo, 
  marker_genes = sub_markers, 
  cluster_key = "clusters_harmony",
  standardize_expression = FALSE
)

# 2. Assign the vector directly to metadata
# This vector is already ordered by cell, so it fits perfectly
data_colo$sub_cell_types_weighted_not_std <- results$annotation_vector

# 3. Check the ranking report to see the scores for each cluster
print(head(results$top5_report))

res_sum_ann <- summarize_annotations(data_colo, "sub_cell_types_weighted_not_std",print_summary = TRUE)
res_sum_ann
# 4. Final Visualization
DimPlot(data_colo, reduction = "umap_harmony", group.by = "sub_cell_types_weighted_not_std", label = TRUE)




results <- get_weighted_annotation(
  seurat_obj = data_colo, 
  marker_genes = sub_markers, 
  cluster_key = "clusters_harmony",
  standardize_expression = TRUE
)

# 2. Assign the vector directly to metadata
# This vector is already ordered by cell, so it fits perfectly
data_colo$sub_cell_types_weighted <- results$annotation_vector

# 3. Check the ranking report to see the scores for each cluster
print(head(results$top5_report))

res_sum_ann <- summarize_annotations(data_colo, "sub_cell_types_weighted",print_summary = TRUE)
res_sum_ann
# 4. Final Visualization
DimPlot(data_colo, reduction = "umap_harmony", group.by = "sub_cell_types_weighted", label = TRUE)







cell_type_markers <- c(
  "Alpi", "Ces2c", "Slc26a2", "Ceacam1", "Aqp8", # Abs. colonocytes
  "Muc2", "Spink4", "Agr2", "Tff3", # Goblet cells
  "Chga", "Chgb", "Tph1", "Cck", "Scgn", "Scg5", # (EECs) Enteroendocrine
  "Dclk1", "Trpm5", "Avil", "Sh2d6", "Plcg2", # Tuft cells
  "Mki67", "Top2a", "Birc5", "Pcna", "Stmn1", # TA cells
  "Lgr5", "Lrig1", "Ascl2", "Slc12a2", "Smoc2", "Kcnq1", "Gpx2",  "Ephb2", "Bmpr1a", "Hopx", "Sox9" #, # Stem cells
  #"Adipoq", "Plin1", "Plin4" # Adipocytes
)

# Create the DotPlot
dot_plot <- DotPlot(object = data_colo, dot.min = 0.05, cols = "RdBu",# c("lightgray", "blue"),# 'RdBu',
                    features = cell_type_markers, scale = T,
                    group.by = "clusters_harmony") +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 13),
    axis.text.y = element_text(size = 15),
    strip.text = element_blank()
  ) 
# Display the plot
dot_plot
ggsave(file.path(output_dir, "Marker_dotplot_harmony_clusters_colonocytes.png"), 
       plot = dot_plot, width = 11, height = 12, dpi = 300, bg = "white")

p0 <- DimPlot(data_colo, reduction = "umap_harmony", group.by = "clusters_harmony", label = TRUE)
p0
ggsave(file.path(output_dir, "UMAP_Cluster_colonocytes_annotate.png"), plot = p0, width = 9, height = 9, dpi = 300)


# Resolution 2.5 annotation
data_colo$sub_cell_types = data_colo$clusters_harmony
data_colo$sub_cell_types = recode_factor(data_colo$sub_cell_types,
                                      '0' = 'Goblet cells',
                                      '1' = 'Abs. colonocytes',
                                      '2' = 'Abs. colonocytes',
                                      '3' = 'Goblet cells',
                                      '4' = 'Stem cells',
                                      '5' = 'Goblet cells',
                                      '6' = 'Abs. colonocytes',
                                      '7' = 'Goblet cells', 
                                      '8' = 'Abs. colonocytes',
                                      '9' = 'Abs. colonocytes',
                                      '10' = 'TA cells',
                                      '11' = 'Abs. colonocytes',
                                      '12' = 'Abs. colonocytes',
                                      '13' = 'TA cells',
                                      '14' = 'TA cells', 
                                      '15' = 'Abs. colonocytes',
                                      '16' = 'Abs. colonocytes',
                                      '17' = 'TA cells',
                                      '18' = 'Stem cells', # or TA cells
                                      '19' = 'Abs. colonocytes',
                                      '20' = 'Goblet cells',
                                      '21' = 'Abs. colonocytes', 
                                      '22' = 'Goblet cells', 
                                      '23' = 'Abs. colonocytes',
                                      '24' = 'Goblet cells',
                                      '25' = 'Goblet cells', 
                                      '26' = 'Goblet cells',
                                      '27' = 'Abs. colonocytes',  
                                      '28' = 'Abs. colonocytes',
                                      '29' = 'TA cells', 
                                      '30' = 'Abs. colonocytes', 
                                      '31' = 'Abs. colonocytes',
                                      '32' = 'Stem cells', # or TA cells
                                      '33' = 'Abs. colonocytes',
                                      '34' = 'Goblet cells',
                                      '35' = 'Abs. colonocytes',
                                      '36' = 'TA cells', 
                                      '37' = 'Abs. colonocytes',
                                      '38' = 'TA cells',
                                      '39' = 'Abs. colonocytes',
                                      '40' = 'EECs',
                                      '41' = 'Goblet cells',
                                      '42' = 'TA cells',
                                      '43' = 'Goblet cells',
                                      '44' = 'Goblet cells',
                                      '45' = 'Tuft cells',
                                      '46' = 'Tuft cells',
                                      '47' = 'EECs',
                                      '48' = 'Goblet cells',
                                      '49' = 'Abs. colonocytes',
                                      '50' = 'Abs. colonocytes',
                                      '51' = 'Goblet cells',
                                      '52' = 'Abs. colonocytes',
                                      '53' = 'Goblet cells',
                                      '54' = 'Goblet cells',
                                      '55' = 'Goblet cells',
                                      '56' = 'Abs. colonocytes',
                                      '57' = 'Goblet cells',
                                      '58' = 'Stem cells',
                                      '59' = 'Goblet cells',
                                      '60' = 'Abs. colonocytes',
                                      '61' = 'TA cells',
                                      '62' = 'Abs. colonocytes'
                                      
)


DimPlot(object = data_colo, reduction = "umap_harmony", group.by = c("sub_cell_types"), label = TRUE, repel = TRUE)


cell_levels <- c("Abs. colonocytes", "Goblet cells", "EECs", 
                 "Tuft cells", "TA cells", "Stem cells")

data_colo$sub_cell_types <- factor(data_colo$sub_cell_types, levels = cell_levels)
data_colo$sub_cell_types_weighted <- factor(data_colo$sub_cell_types_weighted, levels = cell_levels)
data_colo$sub_cell_types_weighted_not_std <- factor(data_colo$sub_cell_types_weighted_not_std, levels = cell_levels)




DimPlot(object = data_colo, reduction = "umap_harmony", 
        group.by = c("sub_cell_types", "sub_cell_types_weighted", "sub_cell_types_weighted_not_std"), 
        label = F, repel = TRUE)



library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl) # Required for Excel export

plot_cell_proportions <- function(seurat_obj, cluster_col, output_prefix, output_dir) {
  
  # 1. Prepare Data
  metadata <- seurat_obj@meta.data
  make_label <- function(p) { ifelse(p > 2.5, paste0(round(p, 1), "%"), "") }
  
  # --- A. Stats by SampleID ---
  df_sample <- metadata %>%
    group_by(SampleID, !!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(SampleID) %>%
    mutate(percentage = (n / sum(n)) * 100)
  
  # --- B. Stats by Genotype_Diet ---
  df_group <- metadata %>%
    group_by(Genotype_Diet, !!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Genotype_Diet) %>%
    mutate(percentage = (n / sum(n)) * 100)
  
  # --- C. Global Stats ---
  df_global <- metadata %>%
    group_by(!!sym(cluster_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(percentage = (n / sum(n)) * 100)
  
  # 2. Save Data to Excel (Multi-sheet)
  excel_path <- file.path(output_dir, paste0(output_prefix, "_Stats.xlsx"))
  
  # Create a named list where names = Sheet Names
  export_list <- list(
    "By_Sample" = df_sample,
    "By_Group"  = df_group,
    "Global"    = df_global
  )
  
  write_xlsx(export_list, path = excel_path)
  message("Excel saved: ", excel_path)
  
  # 3. Generate Plots (Keeping your previous logic)
  p1 <- ggplot(df_sample, aes(x = SampleID, y = percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
    geom_text(aes(label = make_label(percentage)), position = position_stack(vjust = 0.5), size = 4) +
    labs(title = paste("Proportions by Sample (", cluster_col, ")"), y = "Percentage (%)") +
    custom_plot_theme
  
  p2 <- ggplot(df_group, aes(x = Genotype_Diet, y = percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
    geom_text(aes(label = make_label(percentage)), position = position_stack(vjust = 0.5), size = 4) +
    labs(title = paste("Proportions by Group (", cluster_col, ")"), y = "Percentage (%)") +
    custom_plot_theme
  
  p3 <- ggplot(df_global, aes(x = "Global", y = percentage, fill = !!sym(cluster_col))) +
    geom_bar(stat = "identity", color = "white", width = 0.6) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_stack(vjust = 0.5), size = 5) +
    labs(title = paste("Global Distribution (", cluster_col, ")"), y = "Percentage (%)") +
    custom_plot_theme
  
  # 4. Save Plots
  ggsave(file.path(output_dir, paste0(output_prefix, "_cell_pct_by_Sample.png")), p1, width = 12, height = 8, bg = "white")
  ggsave(file.path(output_dir, paste0(output_prefix, "_cell_pct_by_Group.png")),  p2, width = 10, height = 8, bg = "white")
  ggsave(file.path(output_dir, paste0(output_prefix, "_cell_pct_Global.png")),    p3, width = 8, height = 8, bg = "white")
  
  return(list(plots = list(p1, p2, p3), stats = export_list))
}


# Define your labels to test
#labels_to_test <- c("sub_cell_types", "sub_cell_types_weighted", "sub_cell_types_weighted_not_std")
labels_to_test <- c("sub_cell_types")

# Loop through and generate all plots
for (label in labels_to_test) {
  message("Processing proportions for: ", label)
  plot_cell_proportions(
    seurat_obj = data_colo, 
    cluster_col = label, 
    output_prefix = paste0("Colonocytes_", label), 
    output_dir = output_dir
  )
}





cell_type_markers <- c(
  "Alpi", "Ces2c", "Slc26a2", "Ceacam1", "Aqp8", # Abs. colonocytes
  "Muc2", "Spink4", "Agr2", "Tff3", # Goblet cells
  "Chga", "Chgb", "Tph1", "Cck", "Scgn", "Scg5", # (EECs) Enteroendocrine
  "Dclk1", "Trpm5", "Avil", "Sh2d6", "Plcg2", # Tuft cells
  "Mki67", "Top2a", "Birc5", "Pcna", "Stmn1", # TA cells
  "Lgr5", "Lrig1", "Ascl2", "Slc12a2", "Smoc2", "Kcnq1", "Gpx2",  "Ephb2", "Bmpr1a", "Hopx", "Sox9"# Stem cells
)


# UMAP of broad cell types
p1 <- DimPlot(object = data_colo, reduction = "umap_harmony", 
              group.by = "sub_cell_types", 
              label = FALSE, repel = TRUE) +
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16, face = "bold")
  ) + guides(color = guide_legend(override.aes = list(size = 4)))
p1 
ggsave(file.path(output_dir, "UMAP_sub_cell_types.png"), 
       plot = p1, width = 10, height = 7, dpi = 300, bg = "white")


# Create the DotPlot of broad cell type gene markers
dot_plot <- DotPlot(object = data_colo, dot.min = 0.05, cols = "RdBu", #c("white","red"),
                    features = cell_type_markers, scale = T,
                    group.by = "sub_cell_types") + coord_flip() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 13),
    axis.text.y = element_text(size = 15),
    strip.text = element_blank()
  ) 
# Display the plot
dot_plot
ggsave(file.path(output_dir, "Dotplot_markers_sub_cell_types.png"), 
       plot = dot_plot, width = 7, height = 10, dpi = 300, bg = "white")


data_colo$CellType <- data_colo$sub_cell_types
data_colo$seurat_clusters <- data_colo$clusters_harmony



p1 <- DimPlot(
  object = data_colo, 
  reduction = "umap_harmony", 
  group.by = "sub_cell_types", 
  split.by = "Genotype_Diet",  # This "injects" the variable into the plot layers
  label = FALSE, 
  repel = TRUE
) + 
  # 2. Now facet_wrap will recognize the variable
  facet_wrap(~Genotype_Diet, ncol = 2) + 
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16, face = "bold"),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(size = 14, face = "bold") 
  ) + 
  guides(color = guide_legend(override.aes = list(size = 4)))

# Display
p1
# Update dimensions: For a 2x2 grid, a more square aspect ratio works best
ggsave(file.path(output_dir, "UMAP_colonocytes_sub_types_2x2.png"), 
       plot = p1, width = 12, height = 10, dpi = 300, bg = "white")


# Example for a specific marker gene
p1<- generate_gene_comparison_plots(
  seurat_obj = data_colo,
  score_col = "Pparg",
  group_by = "sub_cell_types",
  x_axis = "Genotype_Diet",
  comparisons = list(c("WT_inulin", "KO_inulin"), 
                     c("WT_cellulose", "KO_cellulose"), 
                     c("WT_cellulose", "WT_inulin")),
  plot_type = "violin",
  output_prefix = "pparg_mean_exp_bar_plots_sub_cell_types_",
  plot_title = "Pparg Expression across Groups",
  y_label = "Log-Normalized Expression",
  fig_width = 10, fig_height = 7  
)

p1




# --- SAVE RESULTS ---
message("Saving combined DecontX object...")
saveRDS(data_colo, file.path(output_dir, "data_wu_project_decontX_doubletfinders_colonocytes.rds"))
#saveRDS(data_colo, file.path(output_dir, "data_wu_project_decontX_doubletfinders_colonocytes.rds"), compress = TRUE)



