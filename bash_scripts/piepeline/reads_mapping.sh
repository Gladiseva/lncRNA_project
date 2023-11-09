#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="reads_mapping"
#SBATCH --cpus-per-task=4
#SBATCH --time=00:50:00
#SBATCH --mem-per-cpu=8G

module add UHTS/Aligner/hisat/2.2.1
module add UHTS/Analysis/samtools/1.10

THREADS=$SLURM_CPUS_PER_TASK 

reads=/data/users/lgladiseva/rna_seq/reads
reference_genome=/data/courses/rnaseq_course/lncRNAs/Project1/references/GRCh38.genome.fa
indexed_genome=/data/users/lgladiseva/rna_seq/reference_genome
mapped_reads=/data/users/lgladiseva/rna_seq/mapped_reads

mkdir -p $mapped_reads

for file in "$reads"/*_R1_001*.fastq; do
    # Extract the replicate name from the file name
    replicate=$(basename "$file" | cut -d'_' -f1,2 | sed 's/_L3//' )

    # Construct the file paths for R1 and R2 reads
    r1_file=$(find "$reads" -name "${replicate}_L3_R1_001"*.fastq -type f)
    r2_file=$(find "$reads" -name "${replicate}_L3_R2_001"*.fastq -type f)
    mapped_file="$mapped_reads/${replicate}_aligned_reads.sam"

    # align RNA reads to the genome
    hisat2 -p $THREADS -x "$indexed_genome/GRCh38" -1 "$r1_file" -2 "$r2_file" -S "$mapped_file"
    # sam to bam
    samtools view -@ $THREADS -bS -o "$mapped_reads/${replicate}_aligned_reads.bam" $mapped_file
done



