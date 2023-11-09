#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:55:00
#SBATCH --partition=pall

module add UHTS/Analysis/samtools/1.10

THREADS=$SLURM_CPUS_PER_TASK

samtools view -@ $THREADS -bS -o "/data/users/lgladiseva/rna_seq/mapped_reads/P3_aligned_reads.bam" /data/users/lgladiseva/rna_seq/mapped_reads/P3_L3.sam