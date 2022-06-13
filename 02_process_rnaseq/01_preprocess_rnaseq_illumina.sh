#!/bin/bash
#SBATCH -A snic2018-8-65
#SBATCH -M snowy
#SBATCH -p core -n 4
#SBATCH -t 0-20:00:00
#SBATCH -J trim_galore
#SBATCH -e trim_galore_%A_%a.err            # File to which STDERR will be written
#SBATCH -o trim_galore_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bioinfo-tools TrimGalore/0.4.4

SAMPLEDIR=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/raw_RNAseq_4annotation/190507_A00605_0038_AH7VWFDRXX
WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/process_RNAseq

cd $SAMPLEDIR
FILENAME=`ls -1 *R1*.fastq.gz | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
SAMPLE=$(echo $FILENAME | rev | cut -c 15- | rev | uniq) ; echo $SAMPLE
cd $WORK_D

trim_galore --illumina --paired --retain_unpaired --phred33 --output_dir . --fastqc --gzip --length 36 -q 5 --stringency 1 -e 0.1 ${SAMPLEDIR}/${SAMPLE}1_001.fastq.gz ${SAMPLEDIR}/${SAMPLE}2_001.fastq.gz
#submit array=1-11
