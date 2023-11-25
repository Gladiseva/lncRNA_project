#!/bin/bash

#SBATCH --array=1-12
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=16G

module add UHTS/Quality_control/fastqc/0.11.9

INPUT_DIR="/data/users/lgladiseva/rna_seq/reads"
OUTPUT_DIR="/data/users/lgladiseva/rna_seq/quality_control"

mkdir -p "$OUTPUT_DIR"

# Retrieve the input file for this array task
input_file=$(ls "$INPUT_DIR"/*.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Run FastQC on the selected input file
fastqc -t $SLURM_CPUS_PER_TASK -o "$OUTPUT_DIR" "$input_file"
