#!/bin/bash

#################################################
##### STAR: genome indexing (mouse + human) #####
#################################################

# ------------------------------------------------------------
# -- Run STAR on the latest human and mouse genome versions --
# ------------------------------------------------------------

# Module load
module load STAR/2.7.11a-GCC-13.3.0

# Set timer
start_time=$(date +%s)

# STAR genome indexing (HUMAN)
STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /home/jmartinez/Desktop/toxins_analysis/data/00_reference_genome/indices \
  --genomeFastaFiles /home/jmartinez/Desktop/toxins_analysis/data/00_reference_genome/raw/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile /home/jmartinez/Desktop/toxins_analysis/data/00_reference_genome/raw/gencode.v47.primary_assembly.annotation.gtf \
  --sjdbOverhang 74

# Stop timer
echo "Finished at: $(date)"
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total time taken: $elapsed_time seconds"

