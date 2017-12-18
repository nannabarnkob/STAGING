#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 04.BaseRecalibration.err
#PBS -o 04.BaseRecalibration.log

echo start at `date`

cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load ngs
module load oracle_jdk/1.8.0_144
module load gatk/4.beta.6

runDir=$sample_dir/04.BaseRecalibration
mkdir $runDir
cd $runDir

#java -jar $gatk_launch BaseRecalibrator \
 #   -R $hg38 \
  #  -I $SplitNCigarReads_dir/${sample_name}.split.bam \
   # -knownSites $dbsnp \
    #-knownSites $knownIndels \
    #-L $transcript_intervals \
    #-O $BaseRecalibration_dir/${sample_name}.recal_data.table

java -jar $gatk_launch ApplyBQSR \
    -R $hg38 \
    -I $SplitNCigarReads_dir/${sample_name}.split.bam \
    --bqsr_recal_file $BaseRecalibration_dir/${sample_name}.recal_data.table \
    -O $BaseRecalibration_dir/${sample_name}.recalibrated.bam

echo end at `date`

