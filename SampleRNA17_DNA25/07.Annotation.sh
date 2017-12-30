#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=100gb,walltime=100:00:00
#PBS -e 07.Annotation.err
#PBS -o 07.Annotation.log

echo start at `date`
module load tools
module load anaconda2/4.0.0
module load perl/5.24.0
module load annovar/2016feb01
module load ensembl-tools/89
module load vt/0.5772
module load tabix/1.2.1
module load vcftools/0.1.14

cd $PBS_O_WORKDIR
source Settings.sh

runDir=$sample_dir/07.Annotation
mkdir $runDir
cd $runDir

norm_vcf_suffix=.norm.vcf
annovar_out=.vt.annovar
annovar_protocol=refGene,phastConsElements100way,genomicSuperDups,dbnsfp30a,esp6500siv2_ea,avsnp147,clinvar_20160302,cosmic78,exac03,exac03nontcga,kaviar_20150923,1000g2015aug_eur
vep_out_suffix=.annovar.vep.vcf
vcf_suffix=.vcf.gz

declare -a vcf_names=(${analysis_name}_HC_filtered ${analysis_name}_somaticVariants_filtered);
declare -a paths=($HC_filtered_dir $Mutect2_filtered_dir);

for i in $(seq 0 1)
do
  echo Analyzing ${vcf_names[$i]} ...
  less ${paths[$i]}/${vcf_names[$i]}$vcf_suffix \
 | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
 | vt decompose -s - \
 | vt normalize -r $hg38 - -o $annotation_dir/${vcf_names[i]}$norm_vcf_suffix && \
  perl $annovar $annotation_dir/${vcf_names[$i]}$norm_vcf_suffix $annovar_humandb -buildver hg38 \
    -out $annotation_dir/${vcf_names[$i]}$annovar_out \
    -protocol $annovar_protocol \
    -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput && \
  variant_effect_predictor.pl -i $annotation_dir/${vcf_names[$i]}${annovar_out}.hg38_multianno.vcf \
   --vcf -o $annotation_dir/${vcf_names[$i]}$vep_out_suffix \
   --cache --dir $annotation_cache_dir \
   --sift b --polyphen b --symbol --numbers --biotype --total_length \
   --stats_text --fork 28 \
   --fields Consequence,Distance,Codons,Amino_acids,Gene,SYMBOL,Feature,Feature_type,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE
done
echo All done
echo end at `date


