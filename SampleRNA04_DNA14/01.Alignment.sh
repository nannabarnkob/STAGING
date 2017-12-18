#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 01.Alignment.err
#PBS -o 01.Alignment.log

echo start at `date`

#cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load gcc
module load star/2.5.2b    # alignment of high throughput RNA-seq data - spliced transcripts alignment to a reference)

runDir=$sample_dir/01.Alignment
mkdir $runDir
cd $runDir

STAR --genomeDir $genomeDir \
  --readFilesIn $raw_tumor_R1 $raw_tumor_R2 \
  --readFilesCommand zcat \
  --runThreadN 28 \
  --twopassMode Basic \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix $sample_name"." \
  --outSAMmapqUnique 60

echo end at `date`

