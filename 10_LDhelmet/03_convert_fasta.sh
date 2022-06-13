#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -p core -n 1
#SBATCH -M snowy
#SBATCH -t 0-05:00:00
#SBATCH -J convert
#SBATCH -e convert_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o convert_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#this script converts vcf to the snp format accepted by LDhelmet

ml plink2/2.00-alpha-2-20190429 vcftools/0.1.16 LDhelmet/1.10

WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/03_fasta_format
cd $WORK_D

INTERVAL_LIST=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/Cam_scaffolds.txt
INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'` #31 intervals
SAMPLE_DIR=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/02_shapeit4

#update maf 0.1
vcftools --gzvcf $SAMPLE_DIR/${INTERVAL}.parvulus.phased.vcf.gz --maf 0.1 --out ${INTERVAL}.parvulus --ldhelmet --chr $INTERVAL
