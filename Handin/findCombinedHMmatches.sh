#/bin/sh


# input definition: 
sample_name=$1
# for running manually, use ./findCombinedHMmatches.sh "sample_name" 

module load tools
module load anaconda2/4.0.0


#python /home/projects/dp_00005/data/nanbar/AnnotationScripts/findMatchesKnijnenburg.py ${sample_name}_geneList 
#python /home/projects/dp_00005/data/nanbar/AnnotationScripts/findMatchesNasser.py ${sample_name}_geneList 
python /home/projects/dp_00005/data/nanbar/AnnotationScripts/combineHallmarkTables.py ${sample_name}_HMmatches_Knijnenburg.txt ${sample_name}_HMmatches_Nasser.txt
