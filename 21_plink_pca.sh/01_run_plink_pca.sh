#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 4
#SBATCH -M rackham #snowy was slow for wait time
#SBATCH -t 1-0:00:00
#SBATCH -J plink
#SBATCH -e plink_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o plink_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bioinfo-tools vcftools zlib/1.2.11 plink2

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/plink
cd $WORK_D

INPUT=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/phylogeny/finches_324_shapeit4_phased_rh_GF_prune_autosomes_neutral.vcf
plink2 --vcf $INPUT --pca 4 --out Geospisa_three_species_pruned1 --allow-extra-chr --const-fid --autosome-num 30 --threads 1

###peaks
#
INPUT=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/phylogeny/finches_324_shapeit4_phased_rh_GFautosomes_peaks.vcf
plink2 --vcf $INPUT --pca 4 --out Geospisa_three_species1 --allow-extra-chr --const-fid --autosome-num 30 --threads 4
