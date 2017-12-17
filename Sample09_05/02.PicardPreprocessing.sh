#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 02.PicardPreprocessing.err
#PBS -o 02.PicardPreprocessing.log

echo start at `date`

cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load oracle_jdk/1.8.0_144

runDir=$sample_dir/02.PicardPreprocessing
mkdir $runDir
cd $runDir


java -jar $picard AddOrReplaceReadGroups I=$alignment_dir/${sample_name}Aligned.out.bam O=$PicardProcessing_dir/${sample_name}.rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample 

java -jar $picard MarkDuplicates I=$PicardProcessing_dir/${sample_name}.rg_added_sorted.bam O=$PicardProcessing_dir/${sample_name}.dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 

echo end at `date`

