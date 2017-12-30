#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 06b.FilterMutectCalls.err
#PBS -o 06b.FilterMutectCalls.log

echo start at `date`

cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load ngs
module load oracle_jdk/1.8.0_144
module load gatk/4.beta.6

runDir=$sample_dir/06b.FilterMutectCalls
mkdir $runDir
cd $runDir

java -jar $gatk_launch  FilterMutectCalls \
      -V $Mutect2_dir/${analysis_name}_somaticVariants.vcf.gz \
      -O $runDir/${analysis_name}_somaticVariants_filtered.vcf.gz \
      --dbsnp $dbsnp \
      --normal_artifact_lod 0.0 \
      --tumor_lod 6.0 \
      --base_quality_score_threshold 20 \
      --min_base_quality_score 20 \
      --max_germline_posterior 0.001 \
      --pcr_indel_model HOSTILE \
      --standard_min_confidence_threshold_for_calling 20 \
      --strandArtifactAlleleFraction 0.1 \
      --strandArtifactPosteriorProbability 0.1 \
      --uniqueAltReadCount 5

echo end at `date`

