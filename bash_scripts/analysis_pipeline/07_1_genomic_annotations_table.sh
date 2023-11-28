#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=8G

GENCODE_ALL=/data/courses/rnaseq_course/lncRNAs/Project1/references/gencode.v44.chr_patch_hapl_scaff.annotation.gtf
ALL_MERGED=/data/users/lgladiseva/rna_seq/transcriptome_assembly/ALL_merged.gtf

awk '($11=="transcript_id" && $9=="gene_id") && (/gene_type "protein_coding"/ || /gene_type "lncRNA"/) {print $12,$10,$16,$14}' $GENCODE_ALL | sort -k1,1 | uniq > lncRNA_and_protein_coding.tsv

awk '($11=="transcript_id" && $9=="gene_id" && $13=="gene_name")  {print $12,$10,$14}' $ALL_MERGED | sort -k1,1 | uniq > ALL_genomic_annotations.txt
