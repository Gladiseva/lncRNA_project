#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --partition=pall

module add UHTS/Aligner/hisat/2.2.1

THREADS=$SLURM_CPUS_PER_TASK

READS=/data/users/lgladiseva/rna_seq/reads
INDEXED_GENOME=/data/users/lgladiseva/rna_seq/reference_genome
MAPPED_READS=/data/users/lgladiseva/rna_seq/mapped_reads 

# creates folder if it doesn't exist
mkdir -p /data/users/lgladiseva/rna_seq/mapped_reads 

# align RNA reads to the genome
# -1 and -2 pair reads (forward and reverse),
hisat2 -p $THREADS -x $INDEXED_GENOME/GRCh38 -1 $READS/P3_L3_R1_001_fjv6hlbFgCST.fastq -2 $READS/P3_L3_R2_001_xo7RBLLYYqeu.fastq -S $MAPPED_READS/P3_L3.sam