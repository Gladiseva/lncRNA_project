#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-6

module add UHTS/Analysis/kallisto/0.46.0

THREADS=$SLURM_CPUS_PER_TASK
READS=/data/users/lgladiseva/rna_seq/reads
EXPRESSION=/data/users/lgladiseva/rna_seq/expression_quantification

# Get the current sample ID from the array task ID
current_sample_id=$SLURM_ARRAY_TASK_ID

# Get the corresponding R1 file for the current sample ID
r1_file="$READS/$(ls $READS | grep "_R1_001" | sed -n ${current_sample_id}p)"
# Extract the replicate name from the file name
replicate=$(basename "$r1_file" | cut -d'_' -f1,2 | sed 's/_L3//' )
# Get the corresponding R2 file for the current sample ID
r2_file=$(find "$READS" -name "${replicate}_L3_R2_001"*.fastq -type f)

# Perform Kallisto quantification in paired-end mode
kallisto quant -i $EXPRESSION/all_kallisto_index -o $EXPRESSION/${replicate}_output --rf-stranded --threads $THREADS $r1_file $r2_file
