#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

BED_FILES_LOCATION=/data/users/lgladiseva/rna_seq/inregrative_analysis
BED_FROM_GTF=/data/users/lgladiseva/rna_seq/transcriptome_assembly/ALL_output.bed
STATS_FILE=$BED_FILES_LOCATION/stats_correct_end_annotation.txt

# Count the total number of lines in the original BED file
total_transcripts=$(wc -l < $BED_FROM_GTF)

# Count the number of lines in the overlap files
overlap_5prime=$(wc -l < $BED_FILES_LOCATION/overlap5prime.bed)
overlap_3prime=$(wc -l < $BED_FILES_LOCATION/overlap3prime.bed)

# Calculate the percentages
percentage_5prime=$(bc <<< "scale=2; ($overlap_5prime / $total_transcripts) * 100")
percentage_3prime=$(bc <<< "scale=2; ($overlap_3prime / $total_transcripts) * 100")

# Print the results
echo "Total Transcripts: $total_transcripts" > $STATS_FILE
echo "Transcripts with Correct 5' End Annotation: $overlap_5prime ($percentage_5prime%)" >> $STATS_FILE
echo "Transcripts with Correct 3' End Annotation: $overlap_3prime ($percentage_3prime%)" >> $STATS_FILE


