#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:00:00

module add UHTS/Analysis/BEDTools/2.29.2
module load SequenceAnalysis/GenePrediction/cpat/1.2.4

CPAT_FILES=/data/users/lgladiseva/rna_seq/inregrative_analysis
REFERENCE=/data/courses/rnaseq_course/lncRNAs/Project1/references/GRCh38.genome.fa

bedtools getfasta -fi $REFERENCE -bed $CPAT_FILES/no_overlap53prime.bed -fo $CPAT_FILES/novel_intragenig.fasta

cpat.py -x $CPAT_FILES/Human_Hexamer.tsv -d $CPAT_FILES/Human_logitModel.RData -g $CPAT_FILES/novel_intragenig.fasta -o $CPAT_FILES/novel_intragenic
