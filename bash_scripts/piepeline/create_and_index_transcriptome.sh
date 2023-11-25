#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=04:00:00

module add UHTS/Analysis/kallisto/0.46.0
module add UHTS/Assembler/cufflinks/2.2.1

THREADS=$SLURM_CPUS_PER_TASK
ASSEMBLY=/data/users/lgladiseva/rna_seq/transcriptome_assembly
REFERENCE=/data/courses/rnaseq_course/lncRNAs/Project2/references
EXPRESSION=/data/users/lgladiseva/rna_seq/expression_quantification

mkdir -p $EXPRESSION

#create transcriptome fasta file
gffread -w $EXPRESSION/all_transcriptome.fa -g $REFERENCE/GRCh38.genome.fa $ASSEMBLY/ALL_merged.gtf
# index trnascriptome
kallisto index -i $EXPRESSION/all_kallisto_index $EXPRESSION/all_transcriptome.fa

