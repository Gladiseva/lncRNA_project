#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:50:00

ASSEMBLY=/data/users/lgladiseva/rna_seq/transcriptome_assembly/ALL_merged.gtf
OUTPUT_FILE=/data/users/lgladiseva/rna_seq/transcriptome_assembly/assembly_stats.txt

## Number of:
# Transcripts
transcripts=$(awk '$3=="transcript" {print}' $ASSEMBLY | wc -l)
echo "Number of transcripts: $transcripts" > $OUTPUT_FILE

# Exons
exons=$(awk '$3=="exon" {print}' $ASSEMBLY | wc -l)
echo "Number of exons: $exons" >> $OUTPUT_FILE

# Genes
genes=$(awk '$3=="transcript" {print $10}' $ASSEMBLY | sort -u | wc -l)
echo "Number of genes: $genes" >> $OUTPUT_FILE

## Number of novel:
# Novel Transcripts
novel_transcripts=$(awk '$3=="transcript" && $12 !~ "ENS" {print}' $ASSEMBLY | sort -u | wc -l)
echo "Number of novel transcripts: $novel_transcripts" >> $OUTPUT_FILE

# Novel Exons
novel_exons=$(awk '$3=="exon" && $12 !~ "ENS" {print}' $ASSEMBLY | sort -u | wc -l)
echo "Number of novel exons: $novel_exons" >> $OUTPUT_FILE

## Number of:
# Single Exon Genes
single_exon_genes=$(awk '$3=="exon" {print $10}' $ASSEMBLY |sort | uniq -c | awk '$1 == 1' |  wc -l)
echo "Number of single exon gene : $single_exon_genes" >> $OUTPUT_FILE
# Single Exon Transcripts
single_exon_transcripts=$(awk '$3 == "exon" {print $12}' $ASSEMBLY | sort | uniq -c | awk '$1 == 1' | wc -l)
echo "Number of single exon transcripts: $single_exon_transcripts" >> $OUTPUT_FILE
