#!/bin/bash

#####################################################################
##### STAR: align reads to reference genome (in parallel jobs) #####
#####################################################################

# -----------------------------------------
# -------- Setup the job in SLURM ---------
# -----------------------------------------

#SBATCH --job-name=STAR_alignment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=500G
#SBATCH --partition=bigmem
#SBATCH --output=/users/genomics/jmartinez/STAR_primarycilia_1_log.out
#SBATCH --nodelist=node17

# ------------------------
# ------ Run STAR  -------
# ------------------------

# Load STAR module
module load STAR/2.7.8a-GCC-10.2.0

# Directories
input_dir="/users/genomics/jmartinez/a_primary_cilia_project/01_fastqs/merged"
output_dir="/users/genomics/jmartinez/a_primary_cilia_project/04_STAR"
genome_indices="/users/genomics/jmartinez/data/00_reference_genomes/human/indices"
start_time=$(date +%s)

# Function to process FASTQ files
STAR_align() {
    local FILE="$1"

    if [[ "$FILE" =~ R2_001.fastq.gz$ ]]; then
        return

    elif [[ "$FILE" =~ R1_001.fastq.gz$ ]]; then
        PAIR_1="$FILE"
        PAIR_2="${FILE/R1_001.fastq.gz/R2_001.fastq.gz}"

        if [[ -f "$PAIR_2" ]]; then
            basename_prefix=$(basename "$PAIR_1" | sed 's/R1_001.fastq.gz/_paired_mapped_/')
            STAR \
              --runThreadN 15 \
              --genomeDir "$genome_indices" \
              --readFilesIn "$PAIR_1" "$PAIR_2" \
              --readFilesCommand zcat \
              --outFileNamePrefix "${output_dir}/${basename_prefix}" \
              --outSAMtype BAM SortedByCoordinate
        else
            echo "Warning: Pair file for ${PAIR_1} not found. Skipping..."
        fi
    fi
}

export -f STAR_align

# Set maximum number of parallel jobs
MAX_JOBS=6

# Find all fastq.gz files in input_dir
find "$input_dir" -type f -name "*fastq.gz" | while read -r FILE; do
    STAR_align "$FILE" &
    
    # Enforce job limit dynamically
    while [[ $(jobs -r | wc -l) -ge $MAX_JOBS ]]; do
        wait -n
    done
done

# Wait for all remaining background jobs to finish
wait

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "FASTQ processing completed!"
echo "Total time taken: $elapsed_time seconds"
