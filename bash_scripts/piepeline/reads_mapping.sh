#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --cpus-per-task=8
#SBATCH --time=05:50:00
#SBATCH --mem-per-cpu=8G

module add UHTS/Aligner/hisat/2.2.1
module add UHTS/Analysis/samtools/1.10

THREADS=$SLURM_CPUS_PER_TASK 

READS=/data/users/lgladiseva/rna_seq/reads
INDEXED_GENOME=/data/users/lgladiseva/rna_seq/reference_genome
MAPPED_READS=/data/users/lgladiseva/rna_seq/mapped_reads

mkdir -p $MAPPED_READS

for file in "$READS"/*_R1_001*.fastq; do
    # Extract the replicate name from the file name
    replicate=$(basename "$file" | cut -d'_' -f1,2 | sed 's/_L3//' )

    # Construct the file path for R1 and R2 reads
    r1_file=$(find "$READS" -name "${replicate}_L3_R1_001"*.fastq -type f)
    r2_file=$(find "$READS" -name "${replicate}_L3_R2_001"*.fastq -type f)
    mapped_file="$MAPPED_READS/${replicate}_aligned_reads.sam"

    # align RNA reads to the genome
    hisat2 -p $THREADS -x "$INDEXED_GENOME/GRCh38" -1 "$r1_file" -2 "$r2_file" -S "$mapped_file" --rna-strandness RF > "$MAPPED_READS/${replicate}_hisat2_output.txt" 2>&1
    # sam to bam
    samtools view -@ $THREADS -b -o "$MAPPED_READS/${replicate}_aligned_reads.bam" $mapped_file
done
