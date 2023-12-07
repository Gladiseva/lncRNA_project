#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:00:00

module add UHTS/Analysis/BEDTools/2.29.2

REFERENCES=/data/courses/rnaseq_course/lncRNAs/Project1/references
INTEGRATIVE_ANALYSIS=/data/users/lgladiseva/rna_seq/inregrative_analysis

mkdir -p $OUTPUT_DIR

## -v option stands for "invert," meaning it selects non-overlapping entries.
# intersect for 5' ends
bedtools intersect -v -a $INTEGRATIVE_ANALYSIS/ALL_merged_transcripts.bed -b $REFERENCES/refTSS_v4.1_human_coordinate.hg38.bed > $INTEGRATIVE_ANALYSIS/no_overlap5prime.bed

# intersect for 3' ends
bedtools intersect -v -a $INTEGRATIVE_ANALYSIS/ALL_merged_transcripts.bed -b $REFERENCES/atlas.clusters.2.0.GRCh38.96.bed > $INTEGRATIVE_ANALYSIS/no_overlap3prime.bed

# intesect for 3'end and 5'ends
bedtools intersect -v -a $INTEGRATIVE_ANALYSIS/ALL_merged_transcripts.bed -b $REFERENCES/refTSS_v4.1_human_coordinate.hg38.bed $REFERENCES/atlas.clusters.2.0.GRCh38.96.bed > $INTEGRATIVE_ANALYSIS/no_overlap53prime.bed