#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --partition=pall

module add UHTS/Aligner/hisat/2.2.1

THREADS=$SLURM_CPUS_PER_TASK

READS=/data/users/lgladiseva/rna_seq/reads
INDEXED_GENOME=/data/users/lgladiseva/rna_seq/reference_genome
MAPPED_READS=/data/users/lgladiseva/rna_seq/mapped_reads

mkdir -p /data/users/lgladiseva/rna_seq/mapped_reads

# align RNA reads to the genome
# -1 and -2 pair reads (forward and reverse)
# 2>&1 because hisat2 outputs some important info as errors
# RF -> R1 forward; R2 reverse
hisat2 -p $THREADS -x $INDEXED_GENOME/GRCh38 -1 $READS/3_2_L3_R1_001_DID218YBevN6.fastq -2 $READS/3_2_L3_R2_001_UPhWv8AgN1X1.fastq -S $MAPPED_READS/3_2_L3_new.sam --rna-strandness RF > $MAPPED_READS/3_2_hisat2_output_new.txt 2>&1