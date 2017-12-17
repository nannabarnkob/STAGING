#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 03.SplitNCigarReads.err
#PBS -o 03.SplitNCigarReads.log

echo start at `date`

cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load oracle_jdk/1.8.0_144
module load ngs
module load gatk/4.beta.6

runDir=$sample_dir/03.SplitNCigarReads
mkdir $runDir
cd $runDir

$gatk_launch SplitNCigarReads -R $hg38 -I $PicardPreprocessing_dir/${sample_name}.dedupped.bam -O $SplitNCigarReads_dir/${sample_name}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

echo end at `date`

