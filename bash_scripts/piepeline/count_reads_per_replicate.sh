#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="count_reads_per_sample"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:50:00
#SBATCH --mem=4G
#SBATCH --partition=pall

INPUT_DIR="/data/users/lgladiseva/rna_seq/reads"
OUTPUT_FILE="/data/users/lgladiseva/rna_seq/count_reads_per_sample_output.txt"

declare -A replicate_counts

exec > "$OUTPUT_FILE"
exec 2>&1

for file in "$INPUT_DIR"/*_R1_001*.fastq; do
    # Extract the replicate name from the file name
    replicate=$(basename "$file" | cut -d'_' -f1,2)

    # Count the reads per replicate just from  R1 file
    read_count=$(awk '/^@/{c++}END{print c}' "$file")

    replicate_counts["$replicate"]="$read_count"
done

for replicate in "${!replicate_counts[@]}"; do
    echo "Replicate: $replicate, Read Count: ${replicate_counts["$replicate"]}"
done
