#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 05a.HaplotypeCaller.err
#PBS -o 05a.HaplotypeCaller.log

echo start at `date`

cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load ngs 
module load oracle_jdk/1.8.0_144
module load gatk/4.beta.6

runDir=$sample_dir/05a.HaplotypeCaller
mkdir $runDir
cd $runDir

java -jar $gatk_launch HaplotypeCaller \
  -R $hg38 \
  -O $runDir/${analysis_name}.vcf.gz \
  -I $BaseRecalibration_dir/${sample_name}.recalibrated.bam \
  --max_alternate_alleles 3 \
  -stand_call_conf 20 \
  -dontUseSoftClippedBases \
  -L $transcript_intervals \
  --readFilter OverclippedReadFilter \

echo end at `date`


# -nct 28 removed since it is no longer supported in GATK4. 
# multithreading has been replaced by spark support
