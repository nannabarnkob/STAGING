#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 05b.Mutect2.err
#PBS -o 05b.Mutect2.log

echo start at `date`

cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load ngs
module load oracle_jdk/1.8.0_144
module load gatk/4.beta.6

runDir=$sample_dir/05b.Mutect2
mkdir $runDir
cd $runDir

java -jar $gatk_launch Mutect2 \
    -I $BaseRecalibration_dir/${sample_name}.recalibrated.bam \
    -tumor $sample_name \
    -I $normal_bam \
    -normal Sample_025 \
    -O $Mutect2_dir/${analysis_name}_somaticVariants.vcf.gz \
    -R $hg38
echo end at `date`

