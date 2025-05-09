#!/bin/bash

#############################################
##### Preprocessing and Quality control #####
#############################################

# ----------------------------
# ------ Preprocessing -------
# ----------------------------

base_dir="/home/jmartinez/Desktop/toxins_analysis/data/01_fastqs"
output_dir="$base_dir/merged"
mkdir -p "$output_dir"

# Sample list
sample=(
    "RPE-DMSO-1_S1"
    "RPE-DMSO-2_S2"
    "RPE-DMSO-3_S3"
    "RPE-p8-1_S4"
    "RPE-p8-2_S5"
    "RPE-p8-3_S6"
    "RPE-p15-1_S7"
    "RPE-p15-2_S8"
    "RPE-p15-3_S9"
)
# Loop through each sample and merge files
for sample in "${sample[@]}"; do
    # Find and merge R1 files
    cat $(find "$base_dir" -type f -name "${sample}_L00[1-4]_R1_001.fastq.gz" | sort) > "$output_dir/${sample}_merged_R1_001.fastq.gz"

    # Find and merge R2 files
    cat $(find "$base_dir" -type f -name "${sample}_L00[1-4]_R2_001.fastq.gz" | sort) > "$output_dir/${sample}_merged_R2_001.fastq.gz"
done


# ----------------------------
# ------ Quality Control -----
# ----------------------------

# Load fastqc and multiqc module
module load Miniconda3/20240927
multiqc --version

# Rename input and output directory
base_dir="/home/jmartinez/Desktop/toxins_analysis/data/01_fastqs/merged"
output_dir="/home/jmartinez/Desktop/toxins_analysis/data/02_qc"

# Run FastQC for all files in the directory
fastqc -o "$output_dir" -t 8 "$base_dir"/*.fastq.gz

# Summaries all reports with Multiqc
multiqc $output_dir -o $output_dir

# -----------------------------
# ------ FastP processing -----
# -----------------------------

# Load fastp
module load fastp/0.24.0

# Create output directory for fastp
fastp_output_dir="/home/jmartinez/Desktop/toxins_analysis/data/03_fastp"
mkdir -p "$fastp_output_dir"

# Loop through each sample and process with fastp
for sample in "${sample[@]}"; do
    fastp \
        -i "$output_dir/${sample}_merged_R1_001.fastq.gz" \
        -I "$output_dir/${sample}_merged_R2_001.fastq.gz" \
        -o "$fastp_output_dir/${sample}_filtered_R1.fastq.gz" \
        -O "$fastp_output_dir/${sample}_filtered_R2.fastq.gz" \
        -h "$fastp_output_dir/${sample}_fastp.html" \
        -j "$fastp_output_dir/${sample}_fastp.json" \
        --thread 8
done
