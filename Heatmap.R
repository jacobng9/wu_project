# Install pheatmap if you don't have it already
# install.packages("pheatmap")

library(Seurat)
library(pheatmap)
library(dplyr)

# --- 1. SET PATHS AND LOAD DATA ---
root_path <- "/Users/aayush/RStudio"
output_dir <- file.path(root_path, "seurat_output")
setwd(root_path)

message("Loading Goblet Seurat object...")
data_goblet <- readRDS(file.path(output_dir, "data_wu_project_goblet_subpops_final.rds"))

# --- 2. ORDER THE GROUPS ---
# We want the columns in a logical order on the heatmap: 
# Controls first (Cellulose), then Experimental (Inulin)
data_goblet$Genotype_Diet <- factor(data_goblet$Genotype_Diet, 
                                    levels = c("WT_cellulose", "KO_cellulose", "WT_inulin", "KO_inulin"))
Idents(data_goblet) <- "Genotype_Diet"

# --- 3. DEFINE THE TARGET GENES ---
# Pulling the key hits from your 2-Way ANOVA Enrichment + our sanity check markers
key_genes <- c(
  # The Target / Master Control
  "Pparg", "Muc2", 
  # Starvation & Lipid/Insulin Metabolism
  "SREBF2", "PCSK9", "EIF2AK3", "IDE",
  # Inflammation & NF-kB / Antiviral Defense
  "AIM2", "IFITM3", "OAS2", "OAS3", "ISG15",
  # Barrier Adhesion & Survival
  "BCL2", "RACK1", "DSP", "Ccn3"
)

# Filter the list to ensure all genes are actually present in your current assay
key_genes <- intersect(key_genes, rownames(data_goblet))

# --- 4. CALCULATE AVERAGE EXPRESSION ---
message("Calculating average expression per group...")
# This calculates the mean exponentiated RNA count for each gene per group
avg_exp <- AverageExpression(data_goblet, features = key_genes, return.seurat = FALSE)$RNA

# --- 5. GENERATE AND SAVE THE HEATMAP ---
# We scale by "row" (Z-score) so that high/low expression is relative to each specific gene's baseline.
# Red = Upregulated in that group, Blue = Downregulated in that group

heatmap_filename <- file.path(output_dir, "Goblet_Interaction_Heatmap.png")

# Open a PNG device to save the plot
png(heatmap_filename, width = 6, height = 8, units = "in", res = 300)

pheatmap(
  mat = avg_exp,
  scale = "row",              # Creates the Z-score gradient
  cluster_cols = FALSE,       # Keep our 4 groups in the exact order we defined
  cluster_rows = TRUE,        # Let the genes cluster by similar expression patterns
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "PPARγ x Diet Interaction in Goblet Cells",
  angle_col = 45,             # Tilt the column labels for readability
  fontsize_row = 12,
  fontsize_col = 14
)

# Close the device
dev.off()

message("Heatmap saved to: ", heatmap_filename)