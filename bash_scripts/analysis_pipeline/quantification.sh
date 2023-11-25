#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G


READS=/data/users/lgladiseva/rna_seq/reads

for file in "$READS"/*_R1_001*.fastq; do
    # Extract the replicate name from the file name
    replicate=$(basename "$file" | cut -d'_' -f1,2 | sed 's/_L3//' )
    sbatch /data/users/lgladiseva/rna_seq/scripts_v2/test_per_file/per_file_quantification.sh "$replicate"
done

