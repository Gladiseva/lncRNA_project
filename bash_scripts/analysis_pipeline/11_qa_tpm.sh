#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=8G

# Define an array of output directories
output_dirs=(
    "/data/users/lgladiseva/rna_seq/expression_quantification/3_2_output"
    "/data/users/lgladiseva/rna_seq/expression_quantification/3_4_output"
    "/data/users/lgladiseva/rna_seq/expression_quantification/3_7_output"
    "/data/users/lgladiseva/rna_seq/expression_quantification/P1_output"
    "/data/users/lgladiseva/rna_seq/expression_quantification/P2_output"
    "/data/users/lgladiseva/rna_seq/expression_quantification/P3_output"
)

# Iterate through output directories
for dir in "${output_dirs[@]}"; do
    # Calculate the sum of TPM values
    sum_tpm=$(awk 'NR > 1 { sum += $5 } END { print sum }' "${dir}/abundance.tsv")

    # Print the sum
    echo "Sum of TPM values for ${dir}: $sum_tpm"
done
