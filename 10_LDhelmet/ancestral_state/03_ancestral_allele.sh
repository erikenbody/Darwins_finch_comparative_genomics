#!/bin/bash
#SBATCH -A snic2018-8-65
#SBATCH -p core -n 1
#SBATCH -M snowy
#SBATCH -t 0-03:00:00
#SBATCH -J anc
#SBATCH -e anc_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o anc_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#this runs really fast didnt need to set it up like this, could be run as a for loop

ml plink/1.90b4.9 bcftools

INTERVAL_LIST=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/Cam_scaffolds.txt
INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/01_ancestral_state
cd $WORK_D

echo $INTERVAL $INTERVAL > $INTERVAL.chr
vcftools --bcf barbados.shapeit4.phased_${INTERVAL}.bcf --plink-tped --out $INTERVAL --chrom-map $INTERVAL.chr

#move the Ancestral_allele_tped.pl from the tools dir to this working dir
#then this just puts the path of the tped for this chromosome here
sed "s/outgroup.tped/${INTERVAL}.tped/g" Ancestral_allele_tped.pl > Ancestral_allele_tped_$INTERVAL.pl
perl Ancestral_allele_tped_${INTERVAL}.pl > $INTERVAL.anc
