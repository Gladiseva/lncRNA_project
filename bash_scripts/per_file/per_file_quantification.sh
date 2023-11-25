#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=04:00:00

module add UHTS/Analysis/kallisto/0.46.0

THREADS=$SLURM_CPUS_PER_TASK
READS=/data/users/lgladiseva/rna_seq/reads
ASSEMBLY=/data/users/lgladiseva/rna_seq/transcriptome_assembly/ALL
REFERENCE=/data/courses/rnaseq_course/lncRNAs/Project2/references
EXPRESSION=/data/users/lgladiseva/rna_seq/expression_quantification
REPLICATE=$1

kallisto quant -i $EXPRESSION/all_kallisto_index2 -o $EXPRESSION/${REPLICATE}_output --rf-stranded --threads $THREADS $READS/${REPLICATE}*
