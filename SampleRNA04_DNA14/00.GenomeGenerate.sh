#!/bin/bash
#PBS -W group_list=dp_00005 -A dp_00005
#PBS -l nodes=1:ppn=28,mem=120gb,walltime=15:00:00
#PBS -e 00.GenomeGenerate.err
#PBS -o 00.GenomeGenerate.log

echo start at `date`

cd $PBS_O_WORKDIR
source Settings.sh

module load tools
module load gcc
module load star/2.5.2b    # alignment of high throughput RNA-seq data - spliced transcripts alignment to a reference)

mkdir $genomeDir
cd $genomeDir

STAR  --runMode genomeGenerate \
  --runThreadN 28 \
  --genomeDir $genomeDir \
  --genomeFastaFiles $hg38 \
  --sjdbGTFfile $gtf \
  --sjdbOverhang 100

echo end at `date`

