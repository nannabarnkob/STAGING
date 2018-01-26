#!/bin/bash

module load tools
module load anaconda2/4.0.0

python /home/projects/dp_00005/data/nanbar/AnnotationScripts/combineHallmarkTables.py  HMmatches_BGDistribution.txt HMDistribution_Nasser.txt 
