library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(writexl)
library(ggpubr)

root_path <- "/Users/aayush/RStudio"
# root_path <- "/Users/jacobng/Research/cai/wu_directory/wu_project"
output_dir <- file.path(root_path, "seurat_output")
setwd(root_path)

message("Loading Goblet-specific RDS file...")
data_goblet <- readRDS(file.path(output_dir, "data_wu_project_goblet_subpops_final.rds"))

p1 <- DimPlot(data_goblet, reduction = "umap", group.by = "goblet_subtypes", label = TRUE, repel = TRUE) +
  ggtitle("Goblet Cell Sub-lineages") +
  theme_minimal()

p2 <- DimPlot(data_goblet, reduction = "umap", group.by = "Genotype_Diet", split.by = "Genotype_Diet", ncol = 2) +
  ggtitle("Distribution by Genotype & Diet") +
  theme_minimal()


meta <- data_goblet@meta.data
df_prop <- meta %>%
  group_by(Genotype_Diet, goblet_subtypes) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Genotype_Diet) %>%
  mutate(percentage = (n / sum(n)) * 100)

p3 <- ggplot(df_prop, aes(x = Genotype_Diet, y = percentage, fill = goblet_subtypes)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(title = "Population Proportions", y = "Percentage (%)", x = "") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p4 <- FeaturePlot(data_goblet, features = c("Mki67", "Muc2", "Pparg"), 
                  ncol = 3, cols = c("lightgrey", "red"))

# --- 4. Combine and Save ---
final_plot <- (p1 | p3) / p2 / p4
final_plot

ggsave(file.path(output_dir, "Goblet_Full_Analysis_Plots.png"), 
       plot = final_plot, width = 16, height = 18, dpi = 300, bg = "white")

# Generate a quick summary table of counts
stats_table <- table(data_goblet$goblet_subtypes, data_goblet$Genotype_Diet)
print(stats_table)

# Save as Excel for your records
library(writexl)
write_xlsx(as.data.frame(stats_table), path = file.path(output_dir, "Goblet_Subtype_Counts.xlsx"))
