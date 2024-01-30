# Annotation and characterization of lncRNAs in Lung Cancer

This project analyzes RNA-seq data from non-small cell lung cancer cell populations. The primary goal is to identify potential candidates for cancer treatment among the annotated lncRNA and novel transcripts.

## Introduction

This project objective was to analyze RNA-seq data from biological samples taken from non-small cell lung cancer cell populations. The data for this analysis was used from the study “Tumor Initiation Capacity and Therapy Resistance Are Differential Features of EMT-Related Subpopulations in the NSCLC Cell Line A549”. 

During this project, we analyzed the differentially expressed protein-coding genes, and also annotated and characterized long noncoding RNAs (lncRNAs) and novel genes. Our main goal was to find genes that could serve as potential candidates for treatment from already annotated lncRNA and novel transcripts identified.

## Tools and Versions

Here is a list of the tools used in this project along with their respective versions:

1. **FASTQC:** Version 0.11.9
   - Used for quality control of reads.

2. **HISAT2:** Version 2.2.1
   - Used for mapping reads to the human genome.

3. **SAMTOOLS:** Version 1.10
   - Used for manipulating SAM/BAM files, including sorting and indexing.

4. **StringTie:** Version 1.3.3b
   - Used for reference-guided transcriptome assembly.

5. **CUFFLINKS:** Version 2.2.1
   - Used for creating a transcriptome fasta file.

6. **KALLISTO:** Version 0.46.0
   - Used for transcript quantification.

7. **Sleuth:** Version 0.30.1
   - Used for differential expression analysis.

8. **BEDTools:** Version 2.29.2
   - Used for working with genomic intervals, including finding overlaps.

9. **CPAT:** Version 1.2.4
   - Used for predicting the coding potential of transcripts.

## Usage
Clone repo
Run scripts step-by-step (with sbatch for .sh)
For .sh scripts run on a cluster with sbatch from /data/users/lgladiseva/rna_seq/analysis_pipeline/
R scripts can be run locally.
1. **Step 1: Read quality and statistics:**
   - 01_count_reads_per_replicate.sh
   - 02_quality_control.sh
2. **Step 2: Read mapping:**
   - 03_genome_indexing.sh
   - 04_reads_mapping.sh
   - 05_sort_and_index_bam.sh
3. **Step 3: Transcriptome assembly:**
   - 06_transcriptome_assembly.sh
   - 07_transcriptome_merge.sh
   - 08_assembly_stats.sh
4. **Step 4: Quantification:**
   - 09_create_and_index_transcriptome.sh
   - 10_quantification.sh
   - 11_qa_tpm.sh
5. **Step 5: Differential expression:**
   - 12_count_genes_and_transcripts.R
   - 13_prepare_annotation_table.R
   - 14_sleuth.R
   - 14_2_comparison_with_paper_data.R
6. **Step 6: Integrative analysis:**
   - 15_gtf_to_bed.sh
   - 16_bedtools.sh
   - 17_cpat.sh
   - 18_integrative_analysis_stats.sh
   - 19_construct_table_for_IA.R
7. **Step 7: Prioritization:**
   - 20_prioritization.R

## Where to look for the intermidiate results
cd /data/users/lgladiseva/rna_seq

./expression_quantification
./inregrative_analysis
./mapped_reads
./quality_control
./statistics
./transcriptome_assembly

## Two final CSVs can be found in
/data/users/lgladiseva/rna_seq/RESULTS
