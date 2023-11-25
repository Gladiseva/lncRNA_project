#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=04:00:00

module add UHTS/Aligner/stringtie/1.3.3b

THREADS=$SLURM_CPUS_PER_TASK
REPLICATE=$1
MODE=$2
FILE=$3
BAM_FILES="/data/users/lgladiseva/rna_seq/mapped_reads"
BAM_FILE=${BAM_FILES}/${REPLICATE}_aligned_reads.sorted.bam
ASSEMBLY_DIR="/data/users/lgladiseva/rna_seq/transcriptome_assembly/${MODE}"
REFERENCE="/data/courses/rnaseq_course/lncRNAs/Project2/references"

mkdir -p "$ASSEMBLY_DIR"

# reference-guided transcriptome assembly
# --rf library stranding, -G reference guide annotation
stringtie -o "$ASSEMBLY_DIR/${REPLICATE}_L3.gtf" -p "$THREADS" --rf -G "$REFERENCE/$FILE" "$BAM_FILE"
