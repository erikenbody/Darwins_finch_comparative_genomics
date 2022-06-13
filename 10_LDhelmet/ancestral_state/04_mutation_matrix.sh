#!/bin/bash
#SBATCH -A snic2018-8-65
#SBATCH -p core -n 6
#SBATCH -M snowy
#SBATCH -t 0-03:00:00
#SBATCH -J mut_matrix
#SBATCH -e mut_matrix_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o mut_matrix_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bioinfo-tools bcftools/1.10

WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/01_ancestral_state
cd $WORK_D

cat *anc > outgroup_ancestral.allele

bcftools concat -f parv_vcfs.txt -O v -o parvulus.vcf

#this script found in tools dir
perl Mutation_matrix.pl parvulus.vcf > parvulus_mutation.matrix

#to normalize mutation matrix, see associated excel file
