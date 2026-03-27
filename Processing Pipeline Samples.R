# =============================================================================
# MASTER SCRIPT FOR SINGLE-CELL ANALYSIS (WU IEC PROJECT) 
# =============================================================================

# --- 1. SETUP: LOAD LIBRARIES AND DEFINE PATHS ---
library(Seurat)
library(harmony)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(patchwork)
library(celda)
library(DoubletFinder)
library(writexl)
library(Matrix)

set.seed(123)
root_path <- "Z:/selim_working_dir/2025_wu_iec_project/r_process"
h5_dir <- file.path(root_path, "h5_files")
metadata_file <- file.path(root_path, "Wu_metadata.xlsx")
output_dir <- file.path(root_path, "seurat_output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# =============================================================================


# --- 2. LOAD, MERGE, AND PREPARE DATA ---
message("--- Starting Step 2: Loading, Merging, and Preparing Data ---")

metadata <- read.xlsx(metadata_file)
metadata$SampleID <- as.character(metadata$SampleID)
seurat_objects_list <- list()
for (i in 1:nrow(metadata)) {
  sample_info <- metadata[i, ]
  sample_id <- sample_info$SampleID
  message(paste("\n--- Processing sample:", sample_id, "---"))
  h5_path <- file.path(h5_dir, sample_id, "sample_filtered_feature_bc_matrix.h5")
  if (!file.exists(h5_path)) {
    warning(paste("H5 file not found for sample:", sample_id, ". Skipping."))
    next
  }
  counts_matrix <- Read10X_h5(h5_path)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = sample_id, min.cells = 5)
  for (col_name in colnames(sample_info)) {
    seurat_obj[[col_name]] <- sample_info[[col_name]]
  }
  seurat_obj$batch_id <- paste(seurat_obj$SampleID, seurat_obj$Batch, sep = "_")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^mt-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 500 & percent.mt <= 30)
  seurat_objects_list[[sample_id]] <- seurat_obj
}
rm(seurat_obj, counts_matrix); gc()

if (length(seurat_objects_list) > 0) {
  data <- merge(x = seurat_objects_list[[1]], y = seurat_objects_list[-1], add.cell.ids = names(seurat_objects_list))
  rm(seurat_objects_list); gc()
} else {
  stop("Fatal Error: No data was loaded.")
}

data <- JoinLayers(data)
message("Layers successfully joined into a single counts matrix.")
gc()
# =============================================================================


# --- 3. POST-MERGE QC AND GENE FILTERING ---
message("\n--- Starting Step 3: Post-Merge QC ---")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-|^mt-")
data$plot.group <- "AllCells"

# Create and save the "before" plot
p_before_all <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "plot.group", ncol = 3, pt.size = 0) 
print(p_before_all)
#ggsave(file.path(output_dir, "qc_violin_before_filtering.svg"), plot = p_before_all, width = 10, height = 6)
ggsave(file.path(output_dir, "qc_violin_before_filtering.png"), plot = p_before_all, width = 10, height = 6, dpi = 300)

# Define and apply your QC thresholds
minFeature <- 500
maxFeature <- 14000
minCount <- 1500
maxCount <- 100000
maxMT <- 20
message(paste("Cells before post-merge QC:", ncol(data)))
data <- subset(data, subset =
                 nFeature_RNA >= minFeature & nFeature_RNA <= maxFeature &
                 nCount_RNA >= minCount & nCount_RNA <= maxCount &
                 percent.mt <= maxMT)
message(paste("Cells after post-merge QC:", ncol(data)))

# Create and save the "after" plot
p_after_all <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "plot.group", ncol = 3, pt.size = 0) 
print(p_after_all)
#ggsave(file.path(output_dir, "qc_violin_after_filtering.svg"), plot = p_after_all, width = 10, height = 6)
ggsave(file.path(output_dir, "qc_violin_after_filtering.png"), plot = p_after_all, width = 10, height = 6, dpi = 300)

# Filter genes
genes_to_keep <- rownames(data)[Matrix::rowSums(GetAssayData(data, layer = "counts") > 0) >= 15]
data <- subset(data, features = genes_to_keep)
message(paste("Number of genes remaining:", nrow(data)))
gc()
# =============================================================================






# =============================================================================
# --- 4. DOUBLET DETECTION AND REMOVAL (WITH ROBUST pK SELECTION) ---
# =============================================================================
message("\n--- Starting Step 4: Doublet Detection (Adaptive pK with Fallback) ---")

# Create a subdirectory for the diagnostic pK plots
pk_plot_dir <- file.path(output_dir, "pk_finder_plots2")
if (!dir.exists(pk_plot_dir)) {
  dir.create(pk_plot_dir)
}

data_list <- SplitObject(data, split.by = "SampleID")
results_list <- list()

for (i in seq_along(data_list)) {
  sample_name <- names(data_list)[i]
  message(paste("\n>>> Processing Sample:", sample_name, "(", i, "of", length(data_list), ")"))
  seu_tmp <- data_list[[i]]
  
  # --- Step A: Pre-process Sample ---
  seu_tmp <- NormalizeData(seu_tmp, verbose = FALSE)
  seu_tmp <- FindVariableFeatures(seu_tmp, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu_tmp <- ScaleData(seu_tmp, verbose = FALSE)
  seu_tmp <- RunPCA(seu_tmp, npcs = 50, verbose = FALSE)
  
  # --- Step B: Find Optimal pK with Sanity Checks ---
  message("--- Finding optimal pK ---")
  sweep.res.list <- paramSweep(seu_tmp, PCs = 1:50, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
  
  # --- Step C: Fallback Logic for pK Selection ---
  reasonable_pk_range <- c(0.01, 0.15)
  fallback_pk <- 0.09
  
  initial_optimal_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  final_pk <- initial_optimal_pk
  
  if (final_pk < min(reasonable_pk_range) || final_pk > max(reasonable_pk_range)) {
    message(paste("  [WARNING] Initial optimal pK (", final_pk, ") is outside the reasonable range [",
                  min(reasonable_pk_range), "-", max(reasonable_pk_range), "]. Attempting fallback."))
    
    bcmvn_filtered <- bcmvn[bcmvn$pK > min(reasonable_pk_range) & bcmvn$pK < max(reasonable_pk_range), ]
    
    if (nrow(bcmvn_filtered) > 0 && any(is.finite(bcmvn_filtered$BCmetric))) {
      constrained_pk <- bcmvn_filtered$pK[which.max(bcmvn_filtered$BCmetric)]
      message(paste("  --> Found a new optimal pK within the constrained range:", constrained_pk))
      final_pk <- constrained_pk
    } else {
      message(paste("  --> No suitable peak found in the reasonable range. Reverting to hardcoded fallback pK:", fallback_pk))
      final_pk <- fallback_pk
    }
  } else {
    message(paste("  Optimal pK (", final_pk, ") is within the reasonable range. Proceeding."))
  }
  
  # --- Step D: Create and Save Enhanced Diagnostic Plot ---
  # This section is now AFTER the final_pk is determined.
  
  # Find the BCmetric value corresponding to the final pK
  # Note: If fallback is used, this point might not exist. The plot will still be useful.
  selected_point_data <- bcmvn[bcmvn$pK == final_pk, ]
  
  pk_plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) +
    geom_line(color = "grey60") +
    geom_point(color = "grey60") +
    geom_vline(xintercept = final_pk, linetype = "dashed", color = "red") +
    # Add the red diamond for the selected point ONLY if it exists in the data
    { if(nrow(selected_point_data) > 0)
      geom_point(data = selected_point_data, aes(x = pK, y = BCmetric), color = "red", size = 4, shape = 18)
    } +
    ggtitle(paste0("pK Finder Plot for ", sample_name),
            subtitle = paste0("Final Selected pK = ", final_pk)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(pk_plot_dir, paste0(sample_name, "_pk_plot.png")), plot = pk_plot, width = 7, height = 5, dpi = 150)
  
  # --- Step E: Doublet Estimation ---
  nExp_val <- round(ncol(seu_tmp) * 0.08)
  
  # --- Step F: Run DoubletFinder with the final, validated pK ---
  seu_tmp <- doubletFinder(
    seu_tmp,
    PCs = 1:50,
    pN = 0.25,
    pK = final_pk,
    nExp = nExp_val,
    sct = FALSE
  )
  
  # --- Step G: Save Results ---
  res_col <- grep("^DF.classifications", colnames(seu_tmp@meta.data), value = TRUE)
  res_col <- tail(res_col, 1)
  results_list[[sample_name]] <- seu_tmp@meta.data[, res_col, drop = FALSE]
  
  rm(seu_tmp, sweep.res.list, sweep.stats, bcmvn); gc()
}



# --- 2. Merge Results ---
message("\nMerging results back to main object...")
# 1. Standardize the column names across all dataframes in the list
standardized_list <- lapply(results_list, function(df) {
  colnames(df) <- "Doublet_Status"
  return(df)
})

# 2. Now rbind will work perfectly
all_res <- do.call(rbind, standardized_list)
rownames(all_res) <- sub("^.*?\\.", "", rownames(all_res))

# 3. Map back to the main object
data$Doublet_Status <- all_res[rownames(data@meta.data), "Doublet_Status"]


# Visualize doublets in UMAP
data <- NormalizeData(data)
data <- FindVariableFeatures(data, nfeatures = 5000)
data <- ScaleData(data)
data <- RunPCA(data, npcs = 50, reduction.name = "pca")
gc()

data <- FindNeighbors(data, dims = 1:50, reduction = "pca", k.param = 15, graph.name = "pca_nn")
data <- FindClusters(data, resolution = 1.0, graph.name = "pca_nn", cluster.name = "clusters_none")
data <- RunUMAP(data, dims = 1:50, reduction = "pca", n.neighbors = 15, reduction.name = "umap_none", n.epochs = 500)
gc()
p1 <- DimPlot(data, reduction = "umap_none", group.by = c("Doublet_Status") )
p1 
ggsave(file.path(output_dir, "UMAP_doublet_plot.png"), plot = p1, width = 16, height = 8, dpi = 300)

# Check the results
message("Breakdown of classifications:")
print(table(data$Doublet_Status, useNA = "always"))
# Check cell count before
message(paste("Total cells before filtering:", ncol(data)))
# Logic: Keep if the status is NOT "Doublet" (this keeps Singlets AND NAs)
# We use %in% or is.na for safety
data <- subset(data, subset = Doublet_Status != "Doublet" | is.na(Doublet_Status))
# Check cell count after
message(paste("Total cells after filtering:", ncol(data)))
# Final cleanup
rm(standardized_list, all_res, results_list); gc()



# --- 5. SAVE RAW AND PREP DECONTX ---
message("\n--- Saving Cleaned Raw Object and Preparing DecontX ---")
# Save the raw (singlets-only) object now so it's safe on disk
saveRDS(data, file.path(output_dir, "data_wu_project_post_doubletfinder_raw.rds"))






data<- readRDS(file.path(output_dir, "data_wu_project_post_doubletfinder_raw.rds"))

data <- NormalizeData(data)
data <- FindVariableFeatures(data, nfeatures = 5000)
data <- ScaleData(data)
data <- RunPCA(data, npcs = 50, reduction.name = "pca")
gc()

data <- FindNeighbors(data, dims = 1:50, reduction = "pca", k.param = 15, graph.name = "pca_nn")
data <- FindClusters(data, resolution = 1.0, graph.name = "pca_nn", cluster.name = "clusters_none")
data <- RunUMAP(data, dims = 1:50, reduction = "pca", n.neighbors = 15, reduction.name = "umap_none", n.epochs = 500)
gc()

p1 <- DimPlot(data, reduction = "umap_none", group.by = c("SampleID", "clusters_none" ) )
p1 
ggsave(file.path(output_dir, "UMAP_doublet_removed_no_decontx_or_harmony.png"), plot = p1, width = 16, height = 8, dpi = 300)

p1 <- DimPlot(data, reduction = "umap_none", group.by = c("clusters_none" ) )
p1 
ggsave(file.path(output_dir, "UMAP_doublet_removed_no_decontx_or_harmony_clusters.png"), plot = p1, width = 9, height = 8, dpi = 300)

# =============================================================================
# (The rest of your script for forking into decontX and raw pipelines follows here)
# =============================================================================
# --- 6. DECONTX PIPELINE: HARMONY vs. NON-HARMONY ---
message("\n--- Starting DecontX Pipeline ---")

# 1. Run DecontX
counts_sparse <- GetAssayData(object = data, layer = "counts")
decontx_results <- decontX(x = counts_sparse) # optional z= data$clusters_none, batch=data$SampleID
data[["RNA"]]$counts <- decontx_results$decontXcounts
data$decontX_contamination <- decontx_results$contamination
data[["RNA"]]$data <- NULL 
rm(decontx_results, counts_sparse); gc()

# 2. Refresh Metadata & QC (Essential after DecontX)
data$nCount_RNA <- colSums(data[["RNA"]]$counts)
data$nFeature_RNA <- colSums(data[["RNA"]]$counts > 0)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-|^mt-")
data$plot.group <- "Post-DecontX_Pre-Filter"

# --- QC PLOT: BEFORE FILTERING ---
p_before <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    group.by = "plot.group", ncol = 3, pt.size = 0)
p_before
ggsave(file.path(output_dir, "qc_decontX_BEFORE_filter.png"), plot = p_before, width = 12, height = 6, dpi = 300)

# Apply Subset
data <- subset(data, subset = 
                 nFeature_RNA >= 500 & nFeature_RNA <= 14000 & 
                 nCount_RNA >= 1500 & nCount_RNA <= 100000 & 
                 percent.mt <= 20)

# --- QC PLOT: AFTER FILTERING ---
data$plot.group <- "Post-DecontX_Post-Filter"
p_after <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                   group.by = "plot.group", ncol = 3, pt.size = 0)
p_after
ggsave(file.path(output_dir, "qc_decontX_AFTER_filter.png"), plot = p_after, width = 12, height = 6, dpi = 300)

# 3. Basic Processing
data <- NormalizeData(data)
data <- FindVariableFeatures(data, nfeatures = 5000)
data <- ScaleData(data)
data <- RunPCA(data, npcs = 50, reduction.name = "pca")
gc()

# --- VERSION 1: NON-HARMONY (Standard PCA) ---
message("Processing Non-Harmony Pipeline...")
data <- FindNeighbors(data, dims = 1:50, reduction = "pca", k.param = 15, graph.name = "pca_nn")
data <- FindClusters(data, resolution = 1.0, graph.name = "pca_nn", cluster.name = "clusters_none")
data <- RunUMAP(data, dims = 1:50, reduction = "pca", n.neighbors = 15, reduction.name = "umap_none", n.epochs = 500)
gc()

# --- VERSION 2: HARMONY (Batch Corrected) ---
message("Processing Harmony Pipeline...")
library(harmony)
data <- RunHarmony(data, group.by.vars = "SampleID", reduction = "pca", reduction.save = "harmony")
gc()

data <- FindNeighbors(data, dims = 1:50, reduction = "harmony", k.param = 15, graph.name = "harmony_nn")
data <- FindClusters(data, resolution = 1.0, graph.name = "harmony_nn", cluster.name = "clusters_harmony")
data <- RunUMAP(data, dims = 1:50, reduction = "harmony", n.neighbors = 15, reduction.name = "umap_harmony", n.epochs = 500)
gc()

# --- VISUAL COMPARISON ---
# Side-by-side comparison of Harmony vs. No Harmony
library(ggplot2)

p1 <- DimPlot(data, reduction = "umap_none", group.by = "SampleID") + ggtitle("Standard PCA (No Harmony)")
p2 <- DimPlot(data, reduction = "umap_harmony", group.by = "SampleID") + ggtitle("Harmony Corrected")

# Combine plots and save
p_comp <- p1 + p2 & theme(legend.position = "bottom")
p_comp
ggsave(file.path(output_dir, "UMAP_Harmony_Comparison.png"), plot = p_comp, width = 16, height = 8, dpi = 300)

# To see the actual clusters side-by-side
p3 <- DimPlot(data, reduction = "umap_none", group.by = "clusters_none", label = TRUE) + ggtitle("Clusters: Standard")
p4 <- DimPlot(data, reduction = "umap_harmony", group.by = "clusters_harmony", label = TRUE) + ggtitle("Clusters: Harmony")

p_clusters <- p3 + p4
p_clusters
ggsave(file.path(output_dir, "UMAP_Cluster_Comparison.png"), plot = p_clusters, width = 16, height = 8, dpi = 300)


# --- SAVE RESULTS ---
message("Saving combined DecontX object...")
saveRDS(data, file.path(output_dir, "data_wu_project_decontX_doubletfinders.rds"))
