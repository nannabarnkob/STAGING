#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=100gb,walltime=100:00:00
#PBS -e WDL-Annotation.err
#PBS -o WDL-Annotation.log


# source paths for annotation# 
cd $PBS_O_WORKDIR
source Settings.sh

ln -s /home/projects/dp_00005/data/normal_tumor_workflow/${sample_name}/cromwell-executions/WGS_normal_RNAseq_tumor_SNV_wf/${wf_name}/call-ApplyRecalibrationFilterForSNPs_normal/execution/${sample_name}_normal.recal.final.g.vcf.gz ${dna_name}_filtered.vcf.gz 
ln -s /home/projects/dp_00005/data/normal_tumor_workflow/${sample_name}/cromwell-executions/WGS_normal_RNAseq_tumor_SNV_wf/${wf_name}/call-ApplyRecalibrationFilterForSNPs_normal/execution/${sample_name}_normal.recal.final.g.vcf.gz.tbi ${dna_name}_filtered.vcf.gz.tbi

ln -s /home/projects/dp_00005/data/normal_tumor_workflow/${sample_name}/cromwell-executions/WGS_normal_RNAseq_tumor_SNV_wf/${wf_name}/call-FilterMutectCalls/execution/${sample_name}_MuTect2_filtered.vcf.gz ${sample_name}_MuTect2_filtered.vcf.gz
ln -s /home/projects/dp_00005/data/normal_tumor_workflow/${sample_name}/cromwell-executions/WGS_normal_RNAseq_tumor_SNV_wf/${wf_name}/call-FilterMutectCalls/execution/${sample_name}_MuTect2_filtered.vcf.gz.tbi ${sample_name}_MuTect2_filtered.vcf.gz.tbi

ln -s /home/projects/dp_00005/data/normal_tumor_workflow/${sample_name}/cromwell-executions/WGS_normal_RNAseq_tumor_SNV_wf/${wf_name}/call-VariantFiltration_RNA/execution/${rna_name}_filtered.vcf.gz ${rna_name}_filtered.vcf.gz
ln -s /home/projects/dp_00005/data/normal_tumor_workflow/${sample_name}/cromwell-executions/WGS_normal_RNAseq_tumor_SNV_wf/${wf_name}/call-VariantFiltration_RNA/execution/${rna_name}_filtered.vcf.gz.tbi ${rna_name}_filtered.vcf.gz.tbi

# ANNOTATION PART BEGINGS
echo start at `date`
module load tools
module load anaconda2/4.0.0
module load perl/5.24.0
module load annovar/2017jul16
module load ensembl-tools/89
module load vt/0.5772
module load tabix/1.2.1
module load vcftools/0.1.14

# Alias'for ANNOVAR 
annovar=/services/tools/annovar/2017jul16/table_annovar.pl
annovar_humandb=/home/projects/dp_00005/data/nanbar/ANNOVAR/humandb
reference=/home/databases/gatk-legacy-bundles/b37/human_g1k_v37_decoy.fasta
norm_suffix=.norm.vcf
annovar_protocol=refGene,exac03,dbnsfp30a,esp6500siv2_ea,avsnp147,clinvar_20160302,exac03nontcga,kaviar_20150923
# find all .gz files in the directory and construct array from this (should be 3) 
vcf_files=($(ls *.gz | grep -oP '.*?(?=\.vcf.gz)'))
vcf_suffix=.vcf.gz
annovar_out=.vt.annovar
# For VEP 
vep_out_suffix=.vep.vcf
annotation_cache_dir=/home/databases/variant_effect_predictor/.vep/

for i in $(seq 0 2)
do 
  echo Analyzing ${vcf_files[$i]} ... 
  less ${vcf_files[$i]}$vcf_suffix \
   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
   | vt decompose -s - -o ${vcf_files[$i]}${norm_suffix} && \
  perl $annovar ${vcf_files[$i]}${norm_suffix} $annovar_humandb -buildver hg19 \
    -out ${vcf_files[$i]}${annovar_out} \
    -protocol $annovar_protocol \
    -operation gx,f,f,f,f,f,f,f -nastring . -vcfinput
  echo Done with ANNOVAR protocol for "${vcf_files[$i]}", now making list of gene symbols ... 
  echo Working with "${vcf_files[$i]}".vt.annovar.hg19_multianno.vcf
  grep 'PASS' "${vcf_files[$i]}".vt.annovar.hg19_multianno.vcf | tr ';' '\n' | grep -oP 'Gene.refGene=\K.*' | tr ',' '\n' | sort | uniq > "${vcf_files[$i]}"_geneSymbols 
  sed 's/^\|$/"/g' "${vcf_files[$i]}"_geneSymbols | paste -d, -s > "${vcf_files[$i]}"_geneList \
  echo Lists done, now running variant effect predictor for ${vcf_files[$i]} ... 
  variant_effect_predictor.pl -i ${vcf_files[$i]}${annovar_out}.hg19_multianno.vcf \
   --vcf -o ${vcf_files[$i]}$vep_out_suffix \
   --cache --dir $annotation_cache_dir \
   --sift b --polyphen b --symbol --numbers --biotype --total_length \
   --stats_text --fork 28 \
   --fields Consequence,Distance,Codons,Amino_acids,Gene,SYMBOL,Feature,Feature_type,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE 
  echo Done with "${vcf_files[$i]}"
done 

# CALL PYTHON SCRIPT FOR FINDING MATCHES IN HALLMARKS AND MAKING PLOTS 
python /home/projects/dp_00005/data/nanbar/AnnotationScripts/findMatches.py *_geneList

echo end at `date`
