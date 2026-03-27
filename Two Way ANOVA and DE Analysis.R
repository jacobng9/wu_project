# =============================================================================
# SCRIPT FOR PARALLEL & MEMORY-EFFICIENT TWO-WAY ANOVA ON A LARGE SEURAT OBJECT
# =============================================================================

# --- 1. Setup: Load Libraries and Data ---
# This script is designed to be memory-safe and run in parallel.

# Load required libraries
library(Seurat)
library(dplyr)       # For data manipulation
library(broom)       # To tidy model outputs
library(ggplot2)     # For plotting
library(stringr)     # For string manipulation
library(tidyr)       # For reshaping data

### ADDED FOR PARALLELISM ###
library(future)      # The backend for parallel processing
library(furrr)       # A parallel version of purrr's 'map' functions

# Set seed for reproducibility
set.seed(123)

# Define paths and working directory
root_path <- "Z:/selim_working_dir/2025_wu_iec_project/r_process"
# root_path <- "/mnt/SCDC/Optimus/selim_working_dir/2025_wu_iec_project/r_process"
setwd(root_path)

# --- Explicitly open the large Seurat object file ---
print("Loading the Seurat object. This may take a moment...")
seurat_rds_path <- file.path(root_path, "seurat_output", "data_wu_project_decontX_doubletfinders_colonocytes.rds")
data_colo <- readRDS(seurat_rds_path)
print("Seurat object loaded successfully.")
print(data_colo)


# --- 2. Define Analysis Parameters ---
# --- MODIFY THESE AS NEEDED ---
cell_type_to_analyze <- "Stem cells" # e.g., "Enterocyte", "Goblet cell", "Colonocyte"
cell_type_column <- "sub_cell_types" # Metadata column with cell type annotations
# Factors for the ANOVA model
factor1_name <- "Genotype"
factor2_name <- "Diet"
run_interaction <- TRUE # Test for Genotype * Diet interaction
# Define output paths
output_dir <- file.path(root_path, "seurat_output", "anova_results", cell_type_to_analyze)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_csv_path <- file.path(output_dir, paste0("anova_results_", cell_type_to_analyze, ".csv"))


# --- 3. Verify or Prepare Metadata (Smart Approach) ---
print("Verifying required metadata columns ('Genotype', 'Diet')...")

# Check if the required columns already exist in the metadata
if (all(c(factor1_name, factor2_name) %in% colnames(data_colo@meta.data))) {
  
  print("Found existing 'Genotype' and 'Diet' columns. Verifying format.")
  # The columns exist, so we'll just ensure they are factors, which is required for the ANOVA model.
  # This is a very fast and safe operation.
  data_colo@meta.data <- data_colo@meta.data %>%
    mutate(
      !!factor1_name := as.factor(.data[[factor1_name]]),
      !!factor2_name := as.factor(.data[[factor2_name]])
    )
  # Also verify the interaction column exists, or create it.
  if (!"Genotype_Diet" %in% colnames(data_colo@meta.data)) {
    print("Creating 'Genotype_Diet' column for plotting...")
    data_colo$Genotype_Diet <- paste(data_colo@meta.data[[factor1_name]], data_colo@meta.data[[factor2_name]], sep = "_")
  }
  
} else {
  # If one or more columns are missing, run the original creation code.
  print("One or more required columns not found. Creating them now...")
  data_colo@meta.data <- data_colo@meta.data %>%
    mutate(
      !!factor1_name := case_when(
        str_detect(Batch, "vWT") ~ "WT",
        str_detect(Batch, "vKO") ~ "KO",
        TRUE ~ "Unknown"
      ),
      !!factor2_name := as.factor(Diet)
    )
  # Create the combined 'Genotype_Diet' column for plotting.
  data_colo$Genotype_Diet <- paste(data_colo@meta.data[[factor1_name]], data_colo@meta.data[[factor2_name]], sep = "_")
  print("Metadata columns created.")
}

# This final verification step is still useful to run in either case.
print("Final factor levels for analysis:")
print(table(data_colo@meta.data[[factor1_name]], data_colo@meta.data[[factor2_name]]))



# --- 4. Identify Cells of Interest (NO Subsetting!) ---
# This is the key step for memory efficiency.
print(paste("Identifying cells for analysis in:", cell_type_to_analyze))
Idents(data_colo) <- cell_type_column
cells_to_analyze <- WhichCells(data_colo, idents = cell_type_to_analyze)
metadata_subset <- data_colo@meta.data[cells_to_analyze, c(factor1_name, factor2_name)]
genes_to_test <- rownames(data_colo)
print(paste("Found", length(cells_to_analyze), "cells and will test", length(genes_to_test), "genes."))



# --- 5. Run Two-Way ANOVA in PARALLEL using a CHUNKING STRATEGY with parSapply ---
# This is an alternative implementation using base R's 'parallel' package.
# It demonstrates a different style of parallel programming in R.

# Load the parallel library (usually loaded with R, but good to be explicit)
library(parallel)

# Define the model formula
if (run_interaction) {
  model_formula <- as.formula(paste("expression ~", factor1_name, "*", factor2_name))
} else {
  model_formula <- as.formula(paste("expression ~", factor1_name, "+", factor2_name))
}
print("Using model formula:"); print(model_formula)

### SETUP FOR CHUNKING PARALLELISM with parSapply ###
chunk_size <- 500
gene_chunks <- split(genes_to_test, ceiling(seq_along(genes_to_test) / chunk_size))
num_chunks <- length(gene_chunks)
print(paste("Split", length(genes_to_test), "genes into", num_chunks, "chunks."))

# --- 1. Manually set up the parallel "cluster" ---
num_cores_to_use <- 4
print(paste("Setting up parallel cluster using", num_cores_to_use, "cores..."))
cl <- makeCluster(num_cores_to_use)

# --- 2. Export necessary objects and load libraries on each worker ---
# Each worker needs to know about the model formula and the functions we'll use.
clusterExport(cl, varlist = c("model_formula", "factor1_name", "factor2_name"))
clusterEvalQ(cl, {
  library(broom)
  library(dplyr)
})

all_chunk_results <- list()
print("Starting chunk-based ANOVA processing...")
start_time <- Sys.time()

# --- THE OUTER (SERIAL) LOOP ---
for (i in 1:num_chunks) {
  current_genes <- gene_chunks[[i]]
  print(paste0("--> Processing Chunk ", i, " of ", num_chunks, " (", length(current_genes), " genes)..."))
  
  # A) PRE-FETCH DATA FOR THE CHUNK
  chunk_data_wide <- FetchData(data_colo, vars = c(current_genes, factor1_name, factor2_name), cells = cells_to_analyze)
  
  # B) THE INNER (PARALLEL) LOOP using parLapply (more flexible than parSapply)
  # We pass the cluster object 'cl' to the function.
  chunk_results_list <- parLapply(cl, current_genes, function(gene, chunk_data) {
    gene_df <- data.frame(
      expression = chunk_data[[gene]],
      Genotype   = chunk_data[[factor1_name]],
      Diet       = chunk_data[[factor2_name]]
    )
    aov_model <- aov(model_formula, data = gene_df)
    tidy_results <- broom::tidy(aov_model)
    tidy_results$gene <- gene
    return(tidy_results)
  }, chunk_data = chunk_data_wide) # Pass chunk_data_wide as an extra argument
  
  all_chunk_results[[i]] <- chunk_results_list
} # End of the outer loop

# --- 3. IMPORTANT: Shut down the cluster when you are done ---
stopCluster(cl)



# --- FINALIZING RESULTS ---
anova_results_df <- bind_rows(unlist(all_chunk_results, recursive = FALSE))
end_time <- Sys.time()
print(paste("PARALLEL CHUNKING ANOVA analysis completed in:", round(end_time - start_time, 2), "seconds"))

# The rest of the script remains the same.



# --- 6. Adjust P-values, Format, and SORT the Results Table ---
print("Adjusting p-values, exporting, and sorting results...")

# The name of the interaction term column after pivoting will be "Genotype:Diet.p.adj"
# Let's create this name dynamically so the code is robust.
interaction_term_name <- paste0(factor1_name, ":", factor2_name)
interaction_p_adj_col_name <- paste0(interaction_term_name, ".p.adj")

pvalues_table_sorted <- anova_results_df %>%
  filter(term != "Residuals") %>%
  group_by(term) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  select(gene, term, statistic, p.value, p.adj) %>%
  tidyr::pivot_wider(
    names_from = term,
    values_from = c(statistic, p.value, p.adj),
    names_glue = "{term}.{.value}" # This creates column names like "Genotype:Diet.p.adj"
  ) %>%
  
  ### ADDED FOR SORTING ###
  # Now that we have the wide table, sort it by the adjusted interaction p-value.
  # We use `!!sym()` to use the string variable as a column name inside dplyr.
  arrange(!!sym(interaction_p_adj_col_name))

# Let's look at the top genes with the most significant interaction
print("Top 10 genes with the most significant interaction effect:")
print(head(pvalues_table_sorted, 10))

# Save the SORTED table to the CSV file
write.csv(pvalues_table_sorted, output_csv_path, row.names = FALSE)
print(paste("SORTED results saved to:", output_csv_path))


# --- 7. Plot Results for the TOP Gene (Improved 1x2 Faceted Plot) ---

# Check if there are any results to plot
if (nrow(pvalues_table_sorted) == 0) {
  print("No significant results found to plot.")
} else {
  # Get the top gene from our sorted table
  gene_to_plot <- pvalues_table_sorted$gene[2]
  
  print(paste("Plotting the top gene with the most significant interaction:", gene_to_plot))
  
  # --- 1. Extract and Format P-values for a Clean Subtitle ---
  pvalue_row <- pvalues_table_sorted %>% filter(gene == gene_to_plot)
  
  # Dynamically get the column names for p-values
  genotype_p_col <- paste0(factor1_name, ".p.adj")
  diet_p_col <- paste0(factor2_name, ".p.adj")
  interaction_p_col <- paste0(factor1_name, ":", factor2_name, ".p.adj")
  
  # Create a clean, multi-line subtitle with all relevant ANOVA results
  plot_subtitle <- sprintf(
    "2-Way ANOVA Adj. P-values:\n  Genotype = %.2e\n  Diet = %.2e\n  Interaction = %.2e",
    pvalue_row[[genotype_p_col]],
    pvalue_row[[diet_p_col]],
    pvalue_row[[interaction_p_col]]
  )
  
  # --- 2. Create the Plot using ggplot2 for Full Control ---
  
  # Create the temporary subset object (workaround for your older Seurat version)
  print("Creating temporary subset object for plotting...")
  plotting_subset <- subset(data_colo, cells = cells_to_analyze)
  
  # Fetch the data needed for the plot into a clean data frame
  plot_data <- FetchData(
    plotting_subset,
    vars = c(gene_to_plot, "Genotype", "Diet")
  )
  # Rename the gene column for easier use in ggplot
  colnames(plot_data)[1] <- "Expression"
  
  # Build the plot layer by layer
  violin_plot_final <- ggplot(plot_data, aes(x = Diet, y = Expression, fill = Diet)) +
    # Add violin layer
    geom_violin(trim = TRUE, scale = "width") +
    # Add jittered points for individual cells (optional but good)
    geom_jitter(height = 0, width = 0.1, size = 0.1, alpha = 0.3) +
    
    # This is the key: Facet by Genotype to create two panels (KO and WT)
    facet_wrap(~ Genotype) +
    
    # Set the labels and title
    labs(
      title = paste("Expression of", gene_to_plot, "in", cell_type_to_analyze),
      subtitle = plot_subtitle,
      x = "Diet",
      y = "Log-Normalized Expression"
    ) +
    
    # Use a clean theme
    theme_classic() +
    
    # Customize the theme for readability
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
      strip.text = element_text(face = "bold", size = 14), # "KO" and "WT" labels
      axis.title = element_text(size = 12),
      legend.position = "none" # Hide the legend as it's redundant with the x-axis
    )
  
  # Display the final plot
  print(violin_plot_final)
  
  # Define the new output path for the improved plot
  plot_output_path <- file.path(output_dir, paste0(gene_to_plot, "_", cell_type_to_analyze, "_Final_VlnPlot.png"))
  
  # Save the plot
  ggsave(plot_output_path, plot = violin_plot_final, width = 10, height = 7, dpi = 300)
  
  print(paste("Final, improved plot saved to:", plot_output_path))
  
  # Clean up the temporary object
  rm(plotting_subset, plot_data)
  print("Cleaned up temporary plotting object.")
}
