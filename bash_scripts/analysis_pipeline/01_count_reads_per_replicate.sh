#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=8G

INPUT_DIR="/data/users/lgladiseva/rna_seq/reads"
OUTPUT_FILE="/data/users/lgladiseva/rna_seq/statistics/count_reads_per_sample_output.txt"

declare -A replicate_counts

exec > "$OUTPUT_FILE"
exec 2>&1

# Loop through files per each R1 file
for file in "$INPUT_DIR"/*_R1_001*.fastq; do
    # Extract the replicate name from the file name
    replicate=$(basename "$file" | cut -d'_' -f1,2)

    # Count the reads per replicate from R1 file
    read_count=$(awk '/^@/{c++}END{print c}' "$file")

    replicate_counts["$replicate"]="$read_count"
done

# Loop through replicate names
for replicate in "${!replicate_counts[@]}"; do
    echo "Replicate: $replicate, Read Count: ${replicate_counts["$replicate"]}"
done
