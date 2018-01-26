#/bin/sh
# input defintion 
sample_name=$1
# for running manually, use ./makeHeatmap.sh "sample_name"
module load tools
module load anaconda2/4.4.0

python /home/projects/dp_00005/data/nanbar/AnnotationScripts/makeHeatmap.py ${sample_name}_normal_HMmatches_combined.txt ${sample_name}_tumor_HMmatches_combined.txt ${sample_name}_subtracted_HMmatches_combined.txt ${sample_name}_MuTect_HMmatches_combined.txt
