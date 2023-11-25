#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=05:50:00
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-6 # set the array size according to the number of R1 files

module add UHTS/Aligner/hisat/2.2.1
module add UHTS/Analysis/samtools/1.10

THREADS=$SLURM_CPUS_PER_TASK
READS="/data/users/lgladiseva/rna_seq/reads"
INDEXED_GENOME="/data/users/lgladiseva/rna_seq/reference_genome"
MAPPED_READS="/data/users/lgladiseva/rna_seq/mapped_reads"

mkdir -p "$MAPPED_READS"

# Job array index corresponds to replicate
replicate_index=$SLURM_ARRAY_TASK_ID

# Get the list of R1 files
r1_files=("$READS"/*_R1_001*.fastq)

# Get the R1 file for the current replicate from jobs array
r1_file="${r1_files[$((replicate_index - 1))]}"

# Extract the replicate name from the file name
replicate=$(basename "$r1_file" | cut -d'_' -f1,2 | sed 's/_L3//')

# Construct the file path for R2 reads and the mapped file
r2_file=$(find "$READS" -name "${replicate}_L3_R2_001"*.fastq -type f)
mapped_file="$MAPPED_READS/${replicate}_aligned_reads.sam"

# Align RNA reads to the genome
# -1 and -2 pair reads (forward and reverse)
# 2>&1 because hisat2 outputs some important info as errors
# RF -> R1 forward; R2 reverse
hisat2 -p "$THREADS" -x "$INDEXED_GENOME/GRCh38" -1 "$r1_file" -2 "$r2_file" -S "$mapped_file" --rna-strandness RF > "$MAPPED_READS/${replicate}_hisat2_output.txt" 2>&1

# Convert SAM to BAM
samtools view -@ "$THREADS" -b -o "$MAPPED_READS/${replicate}_aligned_reads.bam" "$mapped_file"
