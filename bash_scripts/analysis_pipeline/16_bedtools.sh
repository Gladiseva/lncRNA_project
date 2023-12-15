#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:00:00

module add UHTS/Analysis/BEDTools/2.29.2

REFERENCES=/data/courses/rnaseq_course/lncRNAs/Project1/references
INTEGRATIVE_ANALYSIS=/data/users/lgladiseva/rna_seq/inregrative_analysis

# find intergenic for novel transcripts
bedtools intersect -v -a $INTEGRATIVE_ANALYSIS/transcripts_novel.bed -b $INTEGRATIVE_ANALYSIS/transcripts_annotated.bed > $INTEGRATIVE_ANALYSIS/novel_intergenic.bed

# 5 prime overlaps
bedtools intersect -wa -s -a $INTEGRATIVE_ANALYSIS/novel5window.bed -b $REFERENCES/refTSS_v4.1_human_coordinate.hg38.bed > $INTEGRATIVE_ANALYSIS/overlap5prime.bed

# 3 prime overlaps
bedtools intersect -wa -s -a $INTEGRATIVE_ANALYSIS/novel3window.bed -b $REFERENCES/atlas.clusters.2.0.GRCh38.96.bed > $INTEGRATIVE_ANALYSIS/overlap3prime.bed