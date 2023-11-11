#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=8G

module add UHTS/Analysis/samtools/1.10

THREADS=$SLURM_CPUS_PER_TASK 

reads=/data/users/lgladiseva/rna_seq/reads
reference_genome=/data/courses/rnaseq_course/lncRNAs/Project1/references/GRCh38.genome.fa
indexed_genome=/data/users/lgladiseva/rna_seq/reference_genome
mapped_reads=/data/users/lgladiseva/rna_seq/mapped_reads

for file in "$reads"/*_R1_001*.fastq; do
    # Extract the replicate name from the file name
    replicate=$(basename "$file" | cut -d'_' -f1,2 | sed 's/_L3//' )
    mapped_file="$mapped_reads/${replicate}_aligned_reads.bam"

    # Sort and index the BAM file
    samtools sort -@ $THREADS -o "/data/users/lgladiseva/rna_seq/mapped_reads/${replicate}_aligned_reads.sorted.bam" $mapped_file
    samtools index -@ $THREADS "/data/users/lgladiseva/rna_seq/mapped_reads/${replicate}_aligned_reads.sorted.bam"

    # Get alignment statistics
    samtools flagstat -@ $THREADS "/data/users/lgladiseva/rna_seq/mapped_reads/${replicate}_aligned_reads.sorted.bam" > "/data/users/lgladiseva/rna_seq/mapped_reads/${replicate}_alignment_stats.txt"
done



