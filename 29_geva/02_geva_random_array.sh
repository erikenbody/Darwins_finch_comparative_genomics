#!/bin/bash
#SBATCH -A snic2021-5-251
#SBATCH -p core -n 8
#SBATCH -t 5-00:00:00
#SBATCH -J geva
#SBATCH -e geva_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o geva_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.14 BEDTools/2.29.2

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/geva
cd $WORK_D

GEVA=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/tools/geva

#geva input files
PEAKS=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/phylogeny/autosomes_for_ful_mag_species_lmm_PEAKS.bed
GZVCF=three_geospiza_phased.vcf.gz
RECOMB=recombination_rate_per_chr_LDHELMET
INTERVAL_LIST=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/lowpass/external_imputation_glimpse/full_ordered_chr.list

chr=`cat $INTERVAL_LIST |  awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

#bedtools intersect -v -a $GZVCF -b $PEAKS -wa -header > three_geospiza_phased_neutral.vcf.gz
GZVCF=three_geospiza_phased_neutral.vcf.gz

#ran this stand alone before
#bcftools view three_geospiza_phased_neutral.vcf -O z -o three_geospiza_phased_neutral.vcf.gz
mkdir -p array_outputs
bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.1 || ALT="*"' -r $chr $GZVCF -O v -o array_outputs/${chr}_three_geospiza_phased_neutral.vcf
awk '{gsub(/^chr1A/,"199"); print}' array_outputs/${chr}_three_geospiza_phased_neutral.vcf | awk '{gsub(/^chr4A/,"499"); print}' | awk '{gsub(/^chr/,""); print}'  > array_outputs/${chr}_three_geospiza_phased_neutral_fixed.vcf #remove "chr" from the name
$GEVA/geva_v1beta --vcf array_outputs/${chr}_three_geospiza_phased_neutral_fixed.vcf --map $RECOMB/${chr}_recomb_rate.txt --out array_outputs/geva_${chr}
awk 'NR%10000==1' array_outputs/geva_${chr}.marker.txt | grep -v "MarkerID"| cut -f 3 -d " " > array_outputs/${chr}_trimmed_positions.txt #take a SNP every 10,000 lines
$GEVA/geva_v1beta -t 8 -i array_outputs/geva_${chr}.bin -o array_outputs/geva_${chr} --positions array_outputs/${chr}_trimmed_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt

rm array_outputs/${chr}_three_geospiza_phased_neutral.vcf
