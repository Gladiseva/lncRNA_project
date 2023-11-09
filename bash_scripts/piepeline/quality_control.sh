#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="quality_control"
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00
#SBATCH --mem=4G

# load fastqc module
module add UHTS/Quality_control/fastqc/0.11.9

INPUT_DIR="/data/users/lgladiseva/rna_seq/reads"
OUTPUT_DIR="/data/users/lgladiseva/rna_seq/quality_control"

mkdir -p $OUTPUT_DIR

for file in $INPUT_DIR/*.fastq;
# default thread count is 4
do
    fastqc -t 2 -o $OUTPUT_DIR "$file";
done