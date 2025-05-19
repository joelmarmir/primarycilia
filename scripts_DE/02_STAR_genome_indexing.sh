#!/bin/bash

#########################################
##### STAR: genome indexing (human) #####
#########################################

# ------------------------------------------------------------
# -- Run STAR on the latest human and mouse genome versions --
# ------------------------------------------------------------

# Module load
module load STAR/2.7.8a-GCC-10.2.0

# Set timer
start_time=$(date +%s)

# STAR genome indexing (HUMAN)
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir /users/genomics/jmartinez/a_primary_cilia_project/00_reference_genome/indices \
  --genomeFastaFiles /users/genomics/jmartinez/a_primary_cilia_project/00_reference_genome/raw/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile /users/genomics/jmartinez/a_primary_cilia_project/00_reference_genome/raw/gencode.v47.primary_assembly.annotation.gtf \
  --sjdbOverhang 74

# Stop timer
echo "Finished at: $(date)"
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total time taken: $elapsed_time seconds"