#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00

module add UHTS/Aligner/stringtie/1.3.3b

THREADS=$SLURM_CPUS_PER_TASK
MODE=$1
FILE=$2
ASSEMBLY_DIR=/data/users/lgladiseva/rna_seq/transcriptome_assembly
REFERENCE=/data/courses/rnaseq_course/lncRNAs/Project2/references

# -p number of threads, --rf library stradness, -G reference annotation
stringtie --merge -o "$ASSEMBLY_DIR/${MODE}_merged.gtf" -p "$THREADS" --rf -G "$REFERENCE/$FILE" "$ASSEMBLY_DIR/$MODE"/*
