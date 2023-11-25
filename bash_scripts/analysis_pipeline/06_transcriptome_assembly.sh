#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=4G

THREADS=$SLURM_CPUS_PER_TASK 

READS=/data/users/lgladiseva/rna_seq/reads

# Define an array of file types and their corresponding GTF files
declare -a file_types=("CHR:gencode.v44.annotation.gtf" "ALL:gencode.v44.chr_patch_hapl_scaff.annotation.gtf" "PRI:gencode.v44.primary_assembly.annotation.gtf")

# Iterate through file types
for type_info in "${file_types[@]}"; do
    IFS=':' read -r -a type_array <<< "$type_info"
    type=${type_array[0]}
    gtf=${type_array[1]}

    # Get every replicate name
    for file in "$READS"/*_R1_001*.fastq; do
        # Extract the replicate name from the file name
        replicate=$(basename "$file" | cut -d'_' -f1,2 | sed 's/_L3//' )

        # Submit the job using the replicate name, type, and GTF file
        sbatch /data/users/lgladiseva/rna_seq/analysis_pipeline/per_file_scripts/per_file_transcriptome_assembly.sh "$replicate" "$type" "$gtf"
    done
done