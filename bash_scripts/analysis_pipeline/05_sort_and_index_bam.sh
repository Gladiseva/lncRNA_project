#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=05:50:00
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-6 # set the array size according to the number of R1 files

module add UHTS/Analysis/samtools/1.10

THREADS=$SLURM_CPUS_PER_TASK 
READS=/data/users/lgladiseva/rna_seq/reads
MAPPED_READS=/data/users/lgladiseva/rna_seq/mapped_reads

# Job array index corresponds to replicate
replicate_index=$SLURM_ARRAY_TASK_ID

# Get the list of R1 files
r1_files=("$READS"/*_R1_001*.fastq)

# Get the R1 file for the current replicate from jobs array
r1_file="${r1_files[$((replicate_index - 1))]}"

# Extract the replicate name from the file name
replicate=$(basename "$r1_file" | cut -d'_' -f1,2 | sed 's/_L3//')

# Construct the file path for replicate
mapped_file="$MAPPED_READS/${replicate}_aligned_reads.bam"
sorted_file="$MAPPED_READS/${replicate}_aligned_reads.sorted.bam"

# Sort and index the BAM file
samtools sort -@ $THREADS -o $sorted_file $mapped_file
samtools index -@ $THREADS $sorted_file

