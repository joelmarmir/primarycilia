################################################
######### Gene Set Enrichment Analysis #########
################################################

# Load libraries
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load data
df_ctrl_vs_p8 <- read.csv("/users/genomics/jmartinez/a_primary_cilia_project/06_fc/ctrl_vs_p8_DESeq2_results1.csv", header = TRUE)
df_ctrl_vs_p15 <- read.csv("/users/genomics/jmartinez/a_primary_cilia_project/06_fc/ctrl_vs_p15_DESeq2_results1.csv", header = TRUE)
df_p8_vs_p15 <- read.csv("/users/genomics/jmartinez/a_primary_cilia_project/06_fc/p8_vs_p15_DESeq2_results1.csv", header = TRUE)

# Select relevant columns and process each dataframe
process_data <- function(df) {
  df %>%
    mutate(score = -log10(padj) * sign(log2FoldChange)) %>%
    arrange(desc(score))
}

df_ctrl_vs_p8 <- process_data(df_ctrl_vs_p8)
df_ctrl_vs_p15 <- process_data(df_ctrl_vs_p15)
df_p8_vs_p15 <- process_data(df_p8_vs_p15)

head(df_ctrl_vs_p8)
tail(df_ctrl_vs_p8)

# Function to prepare gene list
prepare_gene_list <- function(df) {
  geneList <- df$score
  names(geneList) <- df$ENSEMBL_short
  sort(na.omit(geneList), decreasing = TRUE)
}

# Process each comparison
geneList_ctrl_vs_p8 <- prepare_gene_list(df_ctrl_vs_p8)
geneList_ctrl_vs_p15 <- prepare_gene_list(df_ctrl_vs_p15)
geneList_p8_vs_p15 <- prepare_gene_list(df_p8_vs_p15)

head(geneList_ctrl_vs_p8)
head(geneList_ctrl_vs_p15)
head(geneList_p8_vs_p15)
    
# Function to perform GSEA and save results
run_gsea <- function(geneList) {
  GSEA_res <- gseGO(gene = geneList,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENSEMBL",
                   ont = "BP",
                   minGSSize = 15,
                   pvalueCutoff = 0.05,
                   eps = 0)
}

# Run GSEA for all comparisons
GSEA_res_ctrl_vs_p8 <- run_gsea(geneList_ctrl_vs_p8)
GSEA_res_ctrl_vs_p15 <- run_gsea(geneList_ctrl_vs_p15)
GSEA_res_p8_vs_p15 <- run_gsea(geneList_p8_vs_p15)

nrow(GSEA_res_ctrl_vs_p8)
nrow(GSEA_res_ctrl_vs_p15)
nrow(GSEA_res_p8_vs_p15)

head(GSEA_res_ctrl_vs_p8)
head(GSEA_res_ctrl_vs_p15)
head(GSEA_res_p8_vs_p15)


# Visualize results
pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/gsea_ctrl_vs_p8_dotplot.pdf", width = 10, height = 8)
print(dotplot(GSEA_res_ctrl_vs_p8, showCategory = 20, split = ".sign") + ggtitle("GSEA Dotplot - ctrl_vs_p8"))
dev.off()

pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/gsea_ctrl_vs_p15_dotplot.pdf", width = 10, height = 8)
print(dotplot(GSEA_res_ctrl_vs_p15, showCategory = 20, split = ".sign") + ggtitle("GSEA Dotplot - ctrl_vs_p15"))
dev.off()

pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/gsea_p8_vs_p15_dotplot.pdf", width = 10, height = 8)
print(dotplot(GSEA_res_p8_vs_p15, showCategory = 20, split = ".sign") + ggtitle("GSEA Dotplot - p8_vs_p15"))
dev.off()

# Save the results
write.csv(GSEA_res_ctrl_vs_p8, 
          file = "/users/genomics/jmartinez/a_primary_cilia_project/08_ontologies/gsea_BP_ctrl_vs_p8.csv")
write.csv(GSEA_res_ctrl_vs_p15, 
          file = "/users/genomics/jmartinez/a_primary_cilia_project/08_ontologies/gsea_BP_ctrl_vs_p15.csv")
write.csv(GSEA_res_p8_vs_p15, 
          file = "/users/genomics/jmartinez/a_primary_cilia_project/08_ontologies/gsea_BP_p8_vs_p15.csv")

 





"""

sorted_positive_values <- expr_matrix[, "Casey_2018_Ascl1_mesc_DESeq2_results1.csv"] %>%
    sort(decreasing = TRUE)
# Print the name of the column being processed
print(col_name)
# Print the number of positive values
print(length(sorted_positive_values))
# Print the first 10 positive values
print(head(sorted_positive_values, 10))
    
# Calculate gene set enrichment analysis
GSEA_res <- gseGO(gene = sorted_positive_values,
                  OrgDb = "org.Mm.eg.db",
                  keyType = "SYMBOL",
                  ont = "BP")

"""