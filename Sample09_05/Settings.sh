#!/bin/bash

# sample info
sample_name=Sample09_05

# useful paths 
sample_dir=/home/projects/dp_00005/data/nanbar/SampleRuns/Sample09_05


# data files 
raw_tumor_R1=/home/projects/dp_00005/data/RNAseq/Staging04_RNA_S8_R1_0005.fastq.gz
raw_tumor_R2=/home/projects/dp_00005/data/RNAseq/Staging04_RNA_S8_R2_0005.fastq.gz
normal_bam=/170424_XE00510_FCHJMWMALXX_L1_8521703005094.recal.bam
normal_bai=/170424_XE00510_FCHJMWMALXX_L1_8521703005094.recal.bam.bai

# resources
genomeDir=$sample_dir/00.GenomeDir
hg38=/home/projects/dp_00005/data/references/hg38/gatk_chr/GenomeBWAIndex/hg38.fa
gtf=/home/projects/dp_00005/data/references/hg38/hg38.refGene.gtf

# tools 
picard_jar=/services/tools/picard-tools/2.6.0/picard.jar
gatk=/services/tools/gatk/3.6/GenomeAnalysisTK.jar
