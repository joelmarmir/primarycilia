#!/bin/bash

#######################################################################
##### Extract the counts of the mapped reads (in parrallell jobs) #####
#######################################################################

# -----------------------------------------
# -------- Setup the job in SLURM ---------
# -----------------------------------------

#SBATCH --job-name=featureCounts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=500G
#SBATCH --partition=bigmem
#SBATCH --output=/users/genomics/jmartinez/featureCounts_primarycilia_2_log.out
#SBATCH --nodelist=node17

# Load Subread module for featruecounts
module load Subread/2.0.3

# Directories
input_dir="/users/genomics/jmartinez/a_primary_cilia_project/04_STAR"
output_dir="/users/genomics/jmartinez/a_primary_cilia_project/05_counts"
human_gtf="/users/genomics/jmartinez/a_primary_cilia_project/00_reference_genome/raw/gencode.v47.primary_assembly.annotation.gtf"

start_time=$(date +%s)

# Function to process FASTQ files
read_count() {
    local FILE="$1"

    if [[ -f "$FILE" ]]; then
        # Count paired-end BAM files    
        featureCounts \
            -T 15 \
            -p \
            --countReadPairs \
            -B \
            -a "${human_gtf}" \
            -o "${output_dir}/$(basename "$FILE" | sed 's/_mapped_/_counts_/').txt" \
            -s 1 \
            "$FILE"
    # Error handling
    else
        echo "Warning: File not found. Skipping $FILE"
    fi
}
export -f read_count

# Set maximum number of parallel jobs
MAX_JOBS=6

# Find all _1_trimmed.fastq.gz files in input_dir
find "$input_dir" -type f -name "*Aligned.sortedByCoord.out.bam" | while read -r FILE; do
    read_count "$FILE" &
    
    # Enforce job limit dynamically
    while [[ $(jobs -r | wc -l) -ge $MAX_JOBS ]]; do
        wait -n
    done
done

# Wait for all remaining background jobs to finish
wait

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "Counts completed!"
echo "Total time taken: $elapsed_time seconds"

