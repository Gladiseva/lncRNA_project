#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

WORKING_DIR=/data/users/lgladiseva/rna_seq/inregrative_analysis

# sort transcripts
awk '{print $4}' $WORKING_DIR/overlap5prime.bed | sort > $WORKING_DIR/transcripts5prime_sorted.txt
awk '{print $4}' $WORKING_DIR/overlap3prime.bed | sort > $WORKING_DIR/transcripts3prime_sorted.txt

# combine result to see where start and end are correct
join -1 1 -2 1 $WORKING_DIR/transcripts5prime_sorted.txt $WORKING_DIR/transcripts3prime_sorted.txt > $WORKING_DIR/all_overlaps5and3.txt

# count lines in transcripts_novel.bed
total=$(awk '{++n} END{print n}' $WORKING_DIR/transcripts_novel.bed)

# count lines in all_overlaps5and3.txt  
transcripts5and3=$(awk '{++n} END{print n}' $WORKING_DIR/all_overlaps5and3.txt)
transcripts5=$(awk '{++n} END{print n}' $WORKING_DIR/overlap5prime.bed)
transcripts3=$(awk '{++n} END{print n}' $WORKING_DIR/overlap3prime.bed)

# count novel intergenic transcripts
# scale indicates precision
novel_intergenic=$(awk '{++n} END{print n}' $WORKING_DIR/novel_intergenic.bed)
result5=$(echo "scale=3 ; ($novel_intergenic * 100 / $total ) " | bc )
echo "Total of novel transcripts: $total" > $WORKING_DIR/statistics.txt
echo "Number of novel intergenic transcripts: $novel_intergenic" >> $WORKING_DIR/statistics.txt
echo "Percentage of novel transcripts considered as intergenic: $result5%" >> $WORKING_DIR/statistics.txt

# calculate percentage
persentage_5and3=$(echo "scale=3 ; ($transcripts5and3 * 100 / $total ) " | bc )
echo "Percentage of novel transcripts with 5' and 3' good annotations: $persentage_5and3%" >> $WORKING_DIR/statistics.txt

percentage_5=$(echo "scale=3 ; ($transcripts5 * 100 / $total ) " | bc )
echo "Percentage of novel transcripts with 5' good annotations: $percentage_5%" >> $WORKING_DIR/statistics.txt

percentage_3=$(echo "scale=3 ; ($transcripts3 * 100 / $total ) " | bc )
echo "Percentage of novel transcripts with  3' good annotations: $percentage_3%" >> $WORKING_DIR/statistics.txt

# What percent of my novel transcripts are protein coding (0.364 cutoff)
protein_coding=$(awk '{if($5 > 0.364) ++n} END{print n}' $WORKING_DIR/novel_coding_potential.dat)
result=$(echo "scale=3 ; ($protein_coding * 100 / $total ) " | bc )
echo "Percentage of novel transcripts that are protein coding, using 0.364 cutoff: $result%" >> $WORKING_DIR/statistics.txt