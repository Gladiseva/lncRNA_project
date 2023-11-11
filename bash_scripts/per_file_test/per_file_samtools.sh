#!/bin/bash

#SBATCH --CPUS-PER-TASK=4
#SBATCH --MEM-PER-CPU=8G
#SBATCH --TIME=00:55:00
#SBATCH --PARTITION=PALL

module add UHTS/Analysis/samtools/1.10

THREADS=$SLURM_CPUS_PER_TASK

# Input and output file paths
MAPPED_FILE="/data/users/lgladiseva/rna_seq/mapped_reads/3_2_L3_new.sam"
BAM_FILE="/data/users/lgladiseva/rna_seq/mapped_reads/3_2_L3_new.bam"
SORTED_BAM="/data/users/lgladiseva/rna_seq/mapped_reads/3_2_aligned_reads_new.sorted.bam"
ALIGNMENT_STATS="/data/users/lgladiseva/rna_seq/mapped_reads/3_2_new_alignment_stats.txt"

# Convert SAM to BAM
samtools view -@ $THREADS -b -o $BAM_FILE $MAPPED_FILE

# Sort and index the BAM file
samtools sort -@ $THREADS -o $SORTED_BAM $BAM_FILE
samtools index -@ $THREADS $SORTED_BAM

# Get alignment statistics
samtools flagstat -@ $THREADS $SORTED_BAM > $ALIGNMENT_STATS
