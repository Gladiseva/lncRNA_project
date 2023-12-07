#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:00:00

module add UHTS/Analysis/BEDTools/2.29.2

REFERENCES=/data/courses/rnaseq_course/lncRNAs/Project1/references
OUTPUT_DIR=/data/users/lgladiseva/rna_seq/inregrative_analysis
BED_FROM_GTF=/data/users/lgladiseva/rna_seq/transcriptome_assembly/ALL_output.bed

mkdir -p $OUTPUT_DIR

# Bedtools intersect for 5' ends
bedtools intersect -v -a $BED_FROM_GTF -b $REFERENCES/refTSS_v4.1_human_coordinate.hg38.bed > $OUTPUT_DIR/overlap5prime.bed

# Bedtools intersect for 3' ends
bedtools intersect -v -a $BED_FROM_GTF -b $REFERENCES/atlas.clusters.2.0.GRCh38.96.bed > $OUTPUT_DIR/overlap3prime.bed