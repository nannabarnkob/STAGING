#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 06a.VariantFiltration.err
#PBS -o 06a.VariantFiltration.log

echo start at `date`

cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load ngs 
module load oracle_jdk/1.8.0_144
module load gatk/4.beta.6

runDir=$sample_dir/06a.VariantFiltration
mkdir $runDir
cd $runDir

java -jar $gatk_launch VariantFiltration \
      -R $hg38 \
      -V $HC_dir/${analysis_name}.vcf.gz \
      -window 35 \
      -cluster 3 \
      -filterName FS \
      -filter "FS > 30.0" \
      -filterName QD \
      -filter "QD < 2.0" \
      -O $runDir/${analysis_name}_HC_filtered.vcf.gz

echo end at `date`


# -nct 28 removed since it is no longer supported in GATK4. 
# multithreading has been replaced by spark support
