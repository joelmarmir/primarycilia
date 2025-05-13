############################################
##### Differential Expression Analysis #####
############################################

# Load libraries
library(DESeq2)
library(dplyr)
library(purrr)
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
        count_data <- count_table[, c(1, 7, 6)]
        colnames(count_data) <- c("geneid", id, "length")
        
        # Store data
        count_data_list[[id]] <- count_data
    } else {
        print(paste("Count file not found or multiple files found for SRR code:", id))  # Troubleshooting print
    }
}

names(count_data_list)
class(count_data_list)
str(count_data_list)

# Merge all count tables by Geneid
counts <- reduce(count_data_list, full_join, by = c("geneid", "length"))
str(counts)

# 2. Identify sample columns
sample_cols <- setdiff(names(counts), c("geneid", "length"))
str(sample_cols)

# 3. Calculate RPK and TPM
rpk <- counts[sample_cols] / (counts$length / 1e3)
scale_factors <- colSums(rpk) / 1e6
tpm <- sweep(rpk, 2, scale_factors, "/")

head(rpk)
head(tpm)
colSums(tpm)

# 4. Build final TPM data.frame
tpm_df <- bind_cols(
  counts %>% select(geneid, length),
  tpm
)

head(tpm_df)

