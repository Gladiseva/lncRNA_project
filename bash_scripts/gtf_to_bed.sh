#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

module add UHTS/Analysis/bedops/2.4.40

GTF=/data/users/lgladiseva/rna_seq/transcriptome_assembly
BED=/data/users/lgladiseva/rna_seq/inregrative_analysis

#create a directory if it doesn't exist yet
mkdir -p /data/users/lgladiseva/rna_seq/inregrative_analysis

awk '$3=="transcript" {print}' $GTF/ALL_merged.gtf > $GTF/ALL_merged_transcript.gtf

convert2bed -i gtf < $GTF/ALL_merged_transcript.gtf > $BED/ALL_merged_transcripts.bed