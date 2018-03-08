#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=1,mem=4gb,walltime=100:00:00
#PBS -e ex_cromwell_test.err
#PBS -o ex_cromwell_test.log

echo start at `date`

wdl_script=/home/projects/dp_00005/data/normal_tumor_workflow/WGS_normal_RNAseq_tumor_SNV.wdl
wdl_inputs=/home/projects/dp_00005/data/normal_tumor_workflow/Sample_13-016PGF/WGS_normal_RNAseq_tumor_SNV.inputs.json

# load modules
module load moab torque
module load tools
module load oracle_jdk/1.8.0_144
module load gcc/6.2.0
module load R/3.2.5
module load mariadb/10.1.23
module load star/2.5.3a

echo Opening tunnel to padawan ... 
~/mysqltunnel start

# -Dcall-caching.enabled=false
java -Dconfig.file=/home/projects/dp_00005/apps/cromwell.conf -jar /home/projects/dp_00005/apps/src/cromwell/target/scala-2.12/cromwell-30-c4b562f-SNAP.jar run $wdl_script -i $wdl_inputs

echo Closing the tunnel ... 
~/mysqltunnel stop

echo end at `date`

# Considerations: 
# - Memory usage 
# Note on memory: Unless you are using read_string, Cromwell won't ever read the input files into the memory itself, so you don't need to worry about that using up your compute
# Both output paths and job metadata are store in the database. For a large workflow or large number of workflows, this various data stored via an in-memory database can build up, increasing the amount of memory required for the Cromwell process to run.
# https://gatkforums.broadinstitute.org/wdl/discussion/9332/how-much-memory-does-cromwell-need-for-input-or-output-files
# I can confirm that soft-links to previous workflows still work, there is not issue about where the jobs were run (soft-link is created in newest WF for Sample_13-016PGF)
# - We should have a prelimenary step which constructs the workflow from library and choice of call set, something like:  
# mkdir $1 
# cd $1
# $ cat path/to/task_library path/to/paired_wf > final_wf.wdl 
