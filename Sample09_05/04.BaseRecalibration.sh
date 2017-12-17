#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 04.BaseRecalibration.err
#PBS -o 04.BaseRecalibration.log

echo start at `date`

cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load oracle_jdk/1.8.0_144

runDir=$sample_dir/04.BaseRecalibration
mkdir $runDir
cd $runDir

java -jar $gatk_launch \ 
    -T BaseRecalibrator \
    -R $hg38 \
    -I $SplitNCigarReads_dir/${sample_name}.split.bam \
    -knownSites $dbsnp \
    -knownSites $knownIndels \
    -o $BaseRecalibration_dir/${sample_name}.recal_data.table

java -jar $gatk_launch \
    -T PrintReads \
    -R $hg38 \
    -I $SplitNCigarReads_dir/${sample_name}.split.bam \
    -BQSR $BaseRecalibration_dir/${sample_name}.recal_data.table \
    -o $BaseRecalibration_dir/${sample_name}.recalibrated.bam

echo end at `date`

