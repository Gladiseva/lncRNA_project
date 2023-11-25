#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00

module add UHTS/Aligner/hisat/2.2.1

THREADS=$SLURM_CPUS_PER_TASK 
REFERENCE=/data/courses/rnaseq_course/lncRNAs/Project1/references
INDEXED_GENOME=/data/users/lgladiseva/rna_seq/reference_genome

# creates the mappgin folder if it doesn't exist yet
mkdir -p /data/users/lgladiseva/rna_seq/reference_genome

# index the genome
hisat2-build -p $THREADS $REFERENCE/GRCh38.genome.fa $INDEXED_GENOME/GRCh38