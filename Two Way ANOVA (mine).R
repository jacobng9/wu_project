# =============================================================================
# SCRIPT FOR PARALLEL 2-WAY ANOVA (aayush)
# =============================================================================

library(Seurat)
library(dplyr)
library(broom)
library(ggplot2)
library(stringr)
library(tidyr)
library(parallel)

set.seed(123)

# --- 1. macOS PATH CONFIGURATION ---
root_path <- "/Users/aayush/RStudio"
# This points to the specialized Goblet RDS we just saved
seurat_rds_path <- file.path(root_path, "seurat_output", "data_wu_project_goblet_subpops_final.rds")
output_base <- file.path(root_path, "seurat_output", "anova_results")

# Change working directory
setwd(root_path)

print("Loading the Goblet Seurat object...")
data_goblet <- readRDS(seurat_rds_path)
print("Seurat object loaded successfully.")

# --- 2. DEFINE ANALYSIS PARAMETERS ---
# We are now analyzing specific sub-types of Goblet cells
cell_type_to_analyze <- "Canonical_Goblet" 
cell_type_column <- "goblet_subtypes" # The column we just created in the last step

factor1_name <- "Genotype"
factor2_name <- "Diet"
run_interaction <- TRUE 

# Define output directory
output_dir <- file.path(output_base, gsub(" ", "_", cell_type_to_analyze))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_csv_path <- file.path(output_dir, paste0("anova_results_", cell_type_to_analyze, ".csv"))

# --- 3. VERIFY METADATA ---
print("Checking metadata columns...")

if (!"Genotype" %in% colnames(data_goblet@meta.data)) {
  print("Creating 'Genotype' column...")
  data_goblet$Genotype <- ifelse(grepl("vWT", data_goblet$Batch), "WT", 
                                 ifelse(grepl("vKO", data_goblet$Batch), "KO", "Unknown"))
}

data_goblet@meta.data <- data_goblet@meta.data %>%
  mutate(
    !!factor1_name := as.factor(.data[[factor1_name]]),
    !!factor2_name := as.factor(.data[[factor2_name]])
  )

print("Metadata verified. Factor levels:")
print(table(data_goblet@meta.data[[factor1_name]], data_goblet@meta.data[[factor2_name]]))

# --- 4. IDENTIFY CELLS ---
print(paste("Identifying cells for analysis in:", cell_type_to_analyze))
Idents(data_goblet) <- cell_type_column
cells_to_analyze <- WhichCells(data_goblet, idents = cell_type_to_analyze)
genes_to_test <- rownames(data_goblet)

# --- 5. RUN PARALLEL ANOVA ---
if (run_interaction) {
  model_formula <- as.formula(paste("expression ~", factor1_name, "*", factor2_name))
} else {
  model_formula <- as.formula(paste("expression ~", factor1_name, "+", factor2_name))
}

chunk_size <- 500
gene_chunks <- split(genes_to_test, ceiling(seq_along(genes_to_test) / chunk_size))
num_chunks <- length(gene_chunks)

num_cores_to_use <- 4
cl <- makeCluster(num_cores_to_use)
clusterExport(cl, varlist = c("model_formula", "factor1_name", "factor2_name"))
clusterEvalQ(cl, {
  library(broom)
  library(dplyr)
})

all_chunk_results <- list()
start_time <- Sys.time()

for (i in 1:num_chunks) {
  current_genes <- gene_chunks[[i]]
  message(paste0("--> Processing Chunk ", i, " of ", num_chunks))
  
  chunk_data_wide <- FetchData(data_goblet, vars = c(current_genes, factor1_name, factor2_name), cells = cells_to_analyze)
  
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
  }, chunk_data = chunk_data_wide)
  
  all_chunk_results[[i]] <- chunk_results_list
}

stopCluster(cl)

# --- 6. FORMAT & SAVE RESULTS ---
anova_results_df <- bind_rows(unlist(all_chunk_results, recursive = FALSE))
interaction_p_col <- paste0(factor1_name, ":", factor2_name, ".p.adj")

pvalues_table_sorted <- anova_results_df %>%
  filter(term != "Residuals") %>%
  group_by(term) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  select(gene, term, statistic, p.value, p.adj) %>%
  tidyr::pivot_wider(
    names_from = term,
    values_from = c(statistic, p.value, p.adj),
    names_glue = "{term}.{.value}"
  ) %>%
  arrange(!!sym(interaction_p_col))

write.csv(pvalues_table_sorted, output_csv_path, row.names = FALSE)
print(paste("Results saved to:", output_csv_path))

# --- 7. PLOT TOP INTERACTION GENE ---
gene_to_plot <- pvalues_table_sorted$gene[1] # Top hit

plot_data <- FetchData(data_goblet, vars = c(gene_to_plot, factor1_name, factor2_name), cells = cells_to_analyze)
colnames(plot_data)[1] <- "Expression"

vln_plot <- ggplot(plot_data, aes(x = Diet, y = Expression, fill = Diet)) +
  geom_violin(trim = TRUE) +
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.2) +
  facet_wrap(~ Genotype) +
  labs(title = paste(gene_to_plot, "in", cell_type_to_analyze),
       subtitle = "Two-Way ANOVA Interaction Effect",
       y = "Log-Normalized Expression") +
  theme_classic()

ggsave(file.path(output_dir, paste0(gene_to_plot, "_VlnPlot.png")), vln_plot, width = 8, height = 6)