#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

BED_FILES_LOCATION=/data/users/lgladiseva/rna_seq/inregrative_analysis
STATS_FILE=$BED_FILES_LOCATION/stats_correct_end_annotation.txt

# Count the total number of lines in the original BED file
total_transcripts=$(wc -l < $BED_FILES_LOCATION/ALL_merged_transcripts.bed)

# Count the number of lines
no_overlap_5prime=$(wc -l < $BED_FILES_LOCATION/no_overlap5prime.bed)
no_overlap_3prime=$(wc -l < $BED_FILES_LOCATION/no_overlap3prime.bed)
no_overlap_total=$(wc -l < $BED_FILES_LOCATION/no_overlap53prime.bed)

# Calculate the percentages
percentage_5prime_no_overlap=$(bc <<< "scale=2; ($no_overlap_5prime / $total_transcripts) * 100")
percentage_3prime_no_overlap=$(bc <<< "scale=2; ($no_overlap_3prime / $total_transcripts) * 100")
percentage_total_no_overlap=$(bc <<< "scale=2; ($no_overlap_total / $total_transcripts) * 100")

# Print the results
echo "Total transcripts: $total_transcripts" > $STATS_FILE
echo "Transcripts with Correct 5' no overlap End Annotation: $no_overlap_5prime ($percentage_5prime_no_overlap%)" >> $STATS_FILE
echo "Transcripts with Correct 3' no overlap End Annotation: $no_overlap_3prime ($percentage_3prime_no_overlap%)" >> $STATS_FILE
echo "Transcripts with correct start/end no overlap total Annotation: $no_overlap_total ($percentage_total_no_overlap%)" >> $STATS_FILE


