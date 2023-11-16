#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --cpus-per-task=8
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=8G

module add UHTS/Analysis/samtools/1.10

THREADS=$SLURM_CPUS_PER_TASK 

READS=/data/users/lgladiseva/rna_seq/reads
MAPPED_READS=/data/users/lgladiseva/rna_seq/mapped_reads

for file in "$READS"/*_R1_001*.fastq; do
    # Extract the replicate name from the file name
    replicate=$(basename "$file" | cut -d'_' -f1,2 | sed 's/_L3//' )
    mapped_file="$MAPPED_READS/${replicate}_aligned_reads.bam"
    sorted_file="$MAPPED_READS/${replicate}_aligned_reads.sorted.bam"

    # Sort and index the BAM file
    samtools sort -@ $THREADS -o $sorted_file $mapped_file
    samtools index -@ $THREADS $sorted_file
done

