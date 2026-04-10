# =============================================================================
# SCRIPT: 2-WAY ANOVA & ENRICHMENT (ALL GOBLET CELLS, 2 CORES)
# =============================================================================

library(Seurat)
library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(parallel)
library(enrichR)
library(writexl)

set.seed(123)

# --- 1. macOS PATH CONFIGURATION ---
root_path <- "/Users/aayush/RStudio"
seurat_rds_path <- file.path(root_path, "seurat_output", "data_wu_project_goblet_subpops_final.rds")
output_base <- file.path(root_path, "seurat_output", "anova_results")
setwd(root_path)

print("Loading the Goblet Seurat object...")
data_goblet <- readRDS(seurat_rds_path)

# --- 2. DEFINE ANALYSIS PARAMETERS ---
cell_type_to_analyze <- "All_Goblets" 
factor1_name <- "Genotype"
factor2_name <- "Diet"

output_dir <- file.path(output_base, cell_type_to_analyze)
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

# --- 4. IDENTIFY CELLS ---
# We are grabbing ALL cells in this object to contrast the whole Goblet population
cells_to_analyze <- colnames(data_goblet)
genes_to_test <- rownames(data_goblet)
print(paste("Running ANOVA on ALL", length(cells_to_analyze), "Goblet cells across", length(genes_to_test), "genes."))

# --- 5. RUN PARALLEL ANOVA (2 PROCESSORS) ---
model_formula <- as.formula(paste("expression ~", factor1_name, "*", factor2_name))

chunk_size <- 500
gene_chunks <- split(genes_to_test, ceiling(seq_along(genes_to_test) / chunk_size))
num_chunks <- length(gene_chunks)

# Lowered to 2 processors to save RAM!
num_cores_to_use <- 2
cl <- makeCluster(num_cores_to_use)
clusterExport(cl, varlist = c("model_formula", "factor1_name", "factor2_name"))
clusterEvalQ(cl, {
  library(broom)
  library(dplyr)
})

all_chunk_results <- list()

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

# --- 7. PATHWAY ENRICHMENT (ENRICHR) ---
print("Running Pathway Enrichment Analysis on Top Interaction Genes...")

# Grab the top 250 genes driven by the Interaction (Genotype * Diet)
top_genes <- pvalues_table_sorted %>%
  filter(!!sym(interaction_p_col) < 0.05) %>%
  head(250) %>%
  pull(gene)

if(length(top_genes) > 0) {
  dbs <- c("GO_Biological_Process_2023", "KEGG_2019_Mouse", "Reactome_2022")
  setEnrichrSite("Enrichr")
  
  enriched <- enrichr(top_genes, dbs)
  
  # Save to Excel
  write_xlsx(enriched, file.path(output_dir, "Goblet_Interaction_Pathways.xlsx"))
  print("Enrichment pathways saved to Excel.")
  
  # Plot Top 10 GO terms
  p_enrich <- plotEnrich(enriched[[1]], showTerms = 10, numChar = 50, y = "Count", orderBy = "P.value") +
    ggtitle("Top GO Biological Processes (Interaction Effect)") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "Top_Pathways_Barplot.png"), p_enrich, width = 10, height = 6)
  print(p_enrich)
} else {
  print("Not enough significant genes found to run enrichment.")
}