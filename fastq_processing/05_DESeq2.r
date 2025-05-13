############################################
##### Differential Expression Analysis #####
############################################

# Load libraries
library(DESeq2)
library(dplyr)
library(readr)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ggplot2)
library(data.table)

#--------------------------------------
#------------ Functions ---------------
#--------------------------------------

############## Data preparation ##############

# Define file paths
phenodata_path <- "/users/genomics/jmartinez/a_primary_cilia_project/primarycilia/fastq_processing/phenodata.csv"
counts_dir <- "/users/genomics/jmartinez/a_primary_cilia_project/05_counts/strand2"

# Read the phenotype table
phenodata <- read.csv(phenodata_path, header = TRUE, stringsAsFactors = TRUE)

# Ensure replicate and condition is a factor
phenodata$Replicate <- factor(phenodata$Replicate)
phenodata$Condition <- factor(phenodata$Condition)

# Extract SRR codes (assuming they are in a specific column, adjust as needed)
# Let's assume the column is named "SRR_code"
ids <- phenodata$name

# Initialize an empty list to store data
count_data_list <- list()

# Iterate over ids and find corresponding count files
for (id in ids) {
    print(paste("Processing sample:", id))  # Troubleshooting print
    
    count_file <- list.files(counts_dir, pattern = paste0(id, ".*txt$"), full.names = TRUE)
    
    if (length(count_file) == 1) {  # Ensure a single match
        print(paste("Found count file:", count_file))  # Troubleshooting print
        
        count_table <- read.delim(count_file, header = TRUE, stringsAsFactors = FALSE, skip = 1)
        
        # Extract relevant columns: Geneid and counts (assumed to be column 7)
        count_data <- count_table[, c(1, 7)]
        colnames(count_data) <- c("Geneid", id)
        
        # Store data
        count_data_list[[id]] <- count_data
    } else {
        print(paste("Count file not found or multiple files found for SRR code:", id))  # Troubleshooting print
    }
}

# Merge all count tables by Geneid
counts <- Reduce(function(x, y) full_join(x, y, by = "Geneid"), count_data_list)

# Set Geneid as row names
rownames(counts) <- counts$Geneid
counts <- counts[, -1]  # Remove Geneid column

# Generate deseq object with counts and condition information
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = phenodata,
                              design = ~Condition)

############## Data visualization ##############

vsd <- vst(dds, blind = TRUE)

# Export PCA plot to PDF
pdf("PCA_plot.pdf")
plotPCA(vsd, intgroup = c("Condition"))
dev.off()

############## Deseq Analysis ##############

# Filter genes with less than 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds_modified <- dds[keep,]

############## Data visualization ##############

vsd2 <- vst(dds_modified, blind = TRUE)

# Export PCA plot to PDF
pdf("PCA_plot_modified.pdf")
plotPCA(vsd2, intgroup = c("Condition"))
dev.off()

############## Deseq Analysis ##############

# Set uninduced as reference level
dds_modified$Condition <- relevel(dds_modified$Condition, ref = "control")

# Run DESeq and save results
dds_modified <- DESeq(dds_modified)
resultsNames(dds_modified)
res_p8 <- results(dds_modified, name = "Condition_p8_vs_control")
res_p15 <- results(dds_modified, name = "Condition_p15_vs_control")
res_p8_vs_p15 <- results(dds_modified, contrast = c("Condition", "p8", "p15"))


# Extract log2FoldChange and ENSEMBL identifier
res_df_p8 <- as.data.frame(res_p8)
res_df_p15 <- as.data.frame(res_p15)
res_p8_vs_p15 <- as.data.frame(res_p8_vs_p15)

# Remove version number from gene names (to extract geneNames)
res_df_p8$ENSEMBL <- rownames(res_df_p8)
res_df_p15$ENSEMBL <- rownames(res_df_p15)
res_p8_vs_p15$ENSEMBL <- rownames(res_p8_vs_p15)
res_df_p8$ENSEMBL_short <- gsub("\\..*", "",row.names(res_df_p8))
res_df_p15$ENSEMBL_short <- gsub("\\..*", "",row.names(res_df_p15))
res_p8_vs_p15$ENSEMBL_short <- gsub("\\..*", "",row.names(res_p8_vs_p15))

# Map ENSEMBL IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = res_df_p8$ENSEMBL_short, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Add gene symbols to the results data frames
res_df_p8$GeneSymbol <- gene_symbols
res_df_p15$GeneSymbol <- gene_symbols
res_p8_vs_p15$GeneSymbol <- gene_symbols

# Reorder columns to place GeneSymbol first
res_df_p8 <- res_df_p8[, c("ENSEMBL", "ENSEMBL_short", "GeneSymbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "baseMean")]
res_df_p15 <- res_df_p15[, c("ENSEMBL", "ENSEMBL_short", "GeneSymbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "baseMean")]
res_p8_vs_p15 <- res_p8_vs_p15[, c("ENSEMBL", "ENSEMBL_short", "GeneSymbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "baseMean")]

# Remove version number from gene symbols
rownames(res_df_p8) <- NULL
rownames(res_df_p15) <- NULL
rownames(res_p8_vs_p15) <- NULL

# Sort files by log2fc desc
res_df_p8 <- res_df_p8 %>% arrange(padj)
res_df_p15 <- res_df_p15 %>% arrange(padj)
res_p8_vs_p15 <- res_p8_vs_p15 %>% arrange(padj)

write.csv(res_df_p8, file = "/users/genomics/jmartinez/a_primary_cilia_project/06_fc/ctrl_vs_p8_DESeq2_results1.csv", row.names = FALSE)
write.csv(res_df_p15, file = "/users/genomics/jmartinez/a_primary_cilia_project/06_fc/ctrl_vs_p15_DESeq2_results1.csv", row.names = FALSE)
write.csv(res_p8_vs_p15, file = "/users/genomics/jmartinez/a_primary_cilia_project/06_fc/p8_vs_p15_DESeq2_results1.csv", row.names = FALSE)

sum(is.na(res_df_p8$padj))
sum(is.na(res_df_p8$pvalue))
dim(res_df_p8)

head(res_df_p8, 20)


#--------------------------------------
#------- Results visualization --------
#--------------------------------------

############## heatmap ##############   

# Select genes with padj > 0.05 and abs(log2FC) > 0.5
selected_genes_p8 <- res_df_p8 %>% filter(res_df_p8$padj < 0.05, abs(res_df_p8$log2FoldChange) > 0.5)
selected_genes_p15 <- res_df_p15 %>% filter(res_df_p15$padj < 0.05, abs(res_df_p15$log2FoldChange) > 0.5)
selected_genes_p8_vs_p15 <- res_p8_vs_p15 %>% filter(res_p8_vs_p15$padj < 0.05, abs(res_p8_vs_p15$log2FoldChange) > 0.5)

# Check the dimensions of the selected genes
dim(selected_genes_p8)
dim(selected_genes_p15)
dim(selected_genes_p8_vs_p15)

# Merge data frames by ENSEMBL ID
selected_genes <- full_join(selected_genes_p8, selected_genes_p15, by = "ENSEMBL_short")

# Remove duplicates
selected_genes <- distinct(selected_genes)

# Name rows with ENSMBL ID
rownames(selected_genes) <- selected_genes$ENSEMBL_short
selected_genes <- selected_genes[, -1]  # Remove ENSEMBL column

rownames(selected_genes_p8) <- selected_genes_p8$ENSEMBL_short
selected_genes_p8 <- selected_genes_p8[, -1]  # Remove ENSEMBL column

rownames(selected_genes_p15) <- selected_genes_p15$ENSEMBL_short
selected_genes_p15 <- selected_genes_p15[, -1]  # Remove ENSEMBL column

# Check the dimensions of the merged data frame
dim(selected_genes)

# Visualize results in a heatmap
# First normalize countsfrom dds object
hmp_mat <- counts(dds_modified, normalized = TRUE)[(rownames(selected_genes_p15)),]

# Get z-score and name samples
hmp_zmat <- t(apply(hmp_mat, 1, scale))
colnames(hmp_zmat) <- rownames(phenodata)

# Export heatmap to PDF
output_file_heatmap <- "heatmap.pdf"
pdf(output_file_heatmap)
pheatmap(hmp_zmat, cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = phenodata[, c("Condition", "Replicate")])
dev.off()


############## Volcano plot ##############
# Create a function to generate volcano plots with the names of the top 20 DEG
create_volcano_plot <- function(res_df, title) {
    res_df <- res_df %>% mutate(Significant = padj < 0.05 & abs(log2FoldChange) > 0.5)
    
    top_genes <- res_df %>% filter(Significant) %>% arrange(padj) %>% slice(1:10) %>% pull(GeneSymbol)
    
    ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
        geom_point(alpha = 0.5) +
        scale_color_manual(values = c("grey", "red")) +
        theme_minimal() +
        labs(title = title, x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
        theme(legend.position = "none") +
        geom_text(aes(label = ifelse(GeneSymbol %in% top_genes, GeneSymbol, "")), hjust = 1.5, vjust = 1.5)
}

# Generate volcano plots for p8 and p15
pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/volcano_ctrl_vs_p8.pdf")
print(create_volcano_plot(res_df_p8, "Volcano Plot - ctrl vs p8"))
dev.off()

pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/volcano_ctrl_vs_p15.pdf")
print(create_volcano_plot(res_df_p15, "Volcano Plot - ctrl vs p15"))
dev.off()

pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/volcano_p8_vs_p15.pdf")
print(create_volcano_plot(res_p8_vs_p15, "Volcano Plot - p8 vs p15"))
dev.off()

#-------------------------------------------
# --------- Gene ontology -----------------
# ------------------------------------------

# Compute GO
GO_results1_ctrl_vs_p8 <- enrichGO(gene = selected_genes_p8$ENSEMBL_short,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP")

GO_results2_ctrl_vs_p15 <- enrichGO(gene = selected_genes_p15$ENSEMBL_short,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENSEMBL",
                        ont = "BP")

GO_results3_p8_vs_p15 <- enrichGO(gene = selected_genes_p8_vs_p15$ENSEMBL_short,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENSEMBL",
                        ont = "BP")


# Illustrate top 10 ontologies
pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/GO_results1_ctrl_vs_p8_barplot.pdf")
barplot(GO_results1_ctrl_vs_p8, showCategory = 10)
dev.off()

pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/GO_results2_ctrl_vs_p15_barplot.pdf")
barplot(GO_results2_ctrl_vs_p15, showCategory = 10)
dev.off()

pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/GO_results3_p8_vs_p15_barplot.pdf")
barplot(GO_results3_p8_vs_p15, showCategory = 10)
dev.off()


########## Plot genes of interest #############3
gene_info <- read_tsv("/users/genomics/jmartinez/a_primary_cilia_project/99_general/Reiter2017_BigCategories.txt")
head(gene_info)

merge_fc <- function(df, label) {
    df %>%
        select(GeneSymbol, log2FoldChange) %>%
        rename(!!label := log2FoldChange)
}

merged_fc <- gene_info %>%
    left_join(merge_fc(res_df_p8, "log2FC_p8"), by = "GeneSymbol") %>%
    left_join(merge_fc(res_df_p15, "log2FC_p15"), by = "GeneSymbol") %>%
    left_join(merge_fc(res_p8_vs_p15, "log2FC_p8_vs_p15"), by = "GeneSymbol")

head(res_p8_vs_p15)
head(res_df_p8)
head(res_df_p15)



















"""
reiter_genes <- read.delim("/users/genomics/jmartinez/a_primary_cilia_project/99_general/Reiter2017_BigCategories.txt")
genes_of_interest <- unique(reiter_genes$Gene)

# Map symbols to ENSEMBL using the same keytype as your DESeq2 rownames
ensembl_ids <- mapIds(org.Hs.eg.db,
                      keys = genes_of_interest,
                      column = "ENSEMBL",
                      keytype = "SYMBOL",
                      multiVals = "first")

# Remove NAs
ensembl_ids <- na.omit(ensembl_ids)

vsd2 <- vst(dds_modified, blind = TRUE)

# From VST-transformed object (good for visualization)
vsd_mat <- assay(vsd2)
rownames(vsd_mat) <- gsub("\\..*", "", rownames(vsd_mat))
vsd_subset <- vsd_mat[rownames(vsd_mat) %in% ensembl_ids, ]

# Or from normalized counts (raw scale)
norm_counts <- counts(dds_modified, normalized = TRUE)
norm_subset <- norm_counts[rownames(norm_counts) %in% ensembl_ids, ]

annotation <- phenodata[, c("Condition", "Replicate")]
rownames(annotation) <- phenodata$name  # or whatever column matches colnames(vsd_subset)

# Heatmap
pdf("/users/genomics/jmartinez/a_primary_cilia_project/07_plots/heatmap1.pdf")
pheatmap(vsd_subset,
         annotation_col = annotation,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE)
dev.off()

dim(vsd_subset)