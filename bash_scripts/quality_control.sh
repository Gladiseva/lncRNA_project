#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="quality_control"
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00
#SBATCH --mem=4G

# load fastqc module
module add UHTS/Quality_control/fastqc/0.11.9

input_dir="/data/users/lgladiseva/rna_seq/reads"
output_dir="/data/users/lgladiseva/rna_seq/quality_control"

mkdir -p $output_dir

for file in $input_dir/*.fastq;
# default thread count is 4
do
    fastqc -t 2 -o $output_dir "$file";
done