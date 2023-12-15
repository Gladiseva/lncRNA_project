#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=04:00:00

GTF=/data/users/lgladiseva/rna_seq/transcriptome_assembly
INTEGRATIVE_ANALYSIS=/data/users/lgladiseva/rna_seq/inregrative_analysis

#create a directory if it doesn't exist yet
mkdir -p /data/users/lgladiseva/rna_seq/inregrative_analysis

# select only transcript
awk '$3=="transcript" {print}' $GTF/ALL_merged.gtf > $GTF/transcripts_ALL.gtf
# transfrom gtf file in bed file
awk '{if($11=="transcript_id") print $1,$4,$5,$12,$6,$7; else print $1,$4,$5,$12,$6,$7}' $GTF/transcripts_ALL.gtf | sed 's/;//g' | sed 's/"//g' | tr ' ' '\t' > $INTEGRATIVE_ANALYSIS/transcripts.bed

# get annotated transcripts
awk 'match($4,!/MSTR/) {print} ' $INTEGRATIVE_ANALYSIS/transcripts.bed > $INTEGRATIVE_ANALYSIS/transcripts_annotated.bed
#get novel transcripts
awk 'match($4,/MSTR/) {print} ' $INTEGRATIVE_ANALYSIS/transcripts.bed > $INTEGRATIVE_ANALYSIS/transcripts_novel.bed

## window of 100 nt
# 5 prime of novel transcripts
awk '{if($6=="+") print $1"\t"($2-50)"\t"($2+50)"\t"$4"\t"$5"\t"$6 ; else print $1"\t"($3-50)"\t"($3+50)"\t"$4"\t"$5"\t"$6}' $INTEGRATIVE_ANALYSIS/transcripts_novel.bed > $INTEGRATIVE_ANALYSIS/tmp_novel5window.bed
#correction of negative number
awk '{if($3<0 || $2<0) print $1"\t"0"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $INTEGRATIVE_ANALYSIS/tmp_novel5window.bed > $INTEGRATIVE_ANALYSIS/novel5window.bed

# 3 prime of novel transcripts
awk '{if($6=="+") print $1"\t"($3-50)"\t"($3+50)"\t"$4"\t"$5"\t"$6; else print $1"\t"($2-50)"\t"($2+50)"\t"$4"\t"$5"\t"$6}' $INTEGRATIVE_ANALYSIS/transcripts_novel.bed > $INTEGRATIVE_ANALYSIS/tmp_novel3window.bed
#correction of negative number
awk '{if($3<0 || $2<0) print $1"\t"0"\t"$3"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $INTEGRATIVE_ANALYSIS/tmp_novel3window.bed > $INTEGRATIVE_ANALYSIS/novel3window.bed

# Remove temporary files
rm "$INTEGRATIVE_ANALYSIS/tmp_novel5window.bed" "$INTEGRATIVE_ANALYSIS/tmp_novel3window.bed"