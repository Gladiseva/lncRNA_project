#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:50:00

ASSEMBLY=/data/users/lgladiseva/rna_seq/transcriptome_assembly/ALL_merged.gtf
OUTPUT_FILE=/data/users/lgladiseva/rna_seq/transcriptome_assembly/assembly_stats.txt

## Number of:
# Transcripts
n_transcripts=$(awk '$3=="transcript" {print}' $ASSEMBLY | wc -l)
echo "Number of transcripts: $n_transcripts" > $OUTPUT_FILE

# Exons
n_exons=$(awk '$3=="exon" {print}' $ASSEMBLY | wc -l)
echo "Number of exons: $n_exons" >> $OUTPUT_FILE

# Genes
n_genes=$(awk '$3=="transcript" {print $10}' $ASSEMBLY | sort -u | wc -l)
echo "Number of genes: $n_genes" >> $OUTPUT_FILE

## Number of novel:
# Transcripts total - Transcripts with reference
n_novel_transcripts=$(awk '$3=="transcript" && $12 !~ "ENS" {print}' $ASSEMBLY | sort -u | wc -l)
echo "Number of novel transcripts: $n_novel_transcripts" >> $OUTPUT_FILE

# Novel Exons
n_novel_exon=$(awk '$3=="exon" && $12 !~ "ENS" {print}' $ASSEMBLY | sort -u | wc -l)
echo "Number of novel exons: $n_novel_exon" >> $OUTPUT_FILE

## Single exon transcripts/genes
# Single Exon Genes
n_single_exon_gene=$(awk '$3=="exon" {print $10}' $ASSEMBLY |sort | uniq -c | awk '$1 == 1' |  wc -l)
echo "Number of single exon gene : $n_single_exon_gene" >> $OUTPUT_FILE
# Single Exon Transcripts
n_single_exon_transcript=$(awk '$3 == "exon" {print $12}' $ASSEMBLY | sort | uniq -c | awk '$1 == 1' | wc -l)
echo "Number of single exon transcripts: $n_single_exon_transcript" >> $OUTPUT_FILE
