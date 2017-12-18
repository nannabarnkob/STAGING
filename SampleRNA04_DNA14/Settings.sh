#!/bin/bash

# sample info
sample_name=SampleRNA04
DNA_match=DNA14
analysis_name=$sample_name"_"$DNA_match

# useful paths 
sample_dir=/home/projects/dp_00005/data/nanbar/SampleRuns/Sample09_05
genomeDir=$sample_dir/00.GenomeDir                          # Step 00
alignment_dir=$sample_dir/01.Alignment                      # Step 01 
PicardPreprocessing_dir=$sample_dir/02.PicardPreprocessing  # Step 02 
SplitNCigarReads_dir=$sample_dir/03.SplitNCigarReads        # Step 03 
BaseRecalibration_dir=$sample_dir/04.BaseRecalibration      # Step 04 
HC_dir=$sample_dir/05a.HaplotypeCaller                      # Step 05a
Mutect2_dir=$sample_dir/05b.Mutect2                         # Step 05b

# data files 
raw_tumor_R1=/home/projects/dp_00005/data/RNAseq/Staging04_RNA_S8_R1_0005.fastq.gz
raw_tumor_R2=/home/projects/dp_00005/data/RNAseq/Staging04_RNA_S8_R2_0005.fastq.gz
normal_bam=/home/projects/dp_00005/data/petbor/014PGF/02.BSQR/170424_XE00510_FCHJMWMALXX_L1_8521703005094.recal.bam
normal_bai=/home/projects/dp_00005/data/petbor/014PGF/02.BSQR/170424_XE00510_FCHJMWMALXX_L1_8521703005094.recal.bam.bai

# resources
# 0. Several steps 
hg38=/home/projects/dp_00005/data/references/hg38/gatk_chr/GenomeBWAIndex/hg38.fa
transcript_intervals=/home/projects/dp_00005/data/references/hg38/RefGene_transcript.bed # origincally for mutect, but used also for HC

# 1. For Genome Generation (Step 0) 
gtf=/home/projects/dp_00005/data/references/hg38/hg38.refGene.gtf

# 2. For BQSR (Step 04) 
dbsnp=/home/databases/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
knownIndels=/home/projects/dp_00005/data/references/hg38/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf

# 5b. For Mutect2 
gnomad=/home/databases/gnomad_data/vcf/exomes/gnomad.exomes.r2.0.1.sites.vcf.gz

# tools 
picard=/services/tools/picard-tools/2.9.1/picard.jar
gatk38_path=/services/tools/gatk/3.8-0/GenomeAnalysisTK.jar
gatk_launch=/services/tools/gatk/4.beta.6/gatk-package-4.beta.6-local.jar

