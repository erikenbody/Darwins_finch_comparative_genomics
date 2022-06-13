#!/bin/bash
#SBATCH -A snic2021-5-251
#SBATCH -p core -n 5
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

#bedtools intersect -v -a $GZVCF -b $PEAKS -wa -header > three_geospiza_phased_neutral.vcf.gz
GZVCF=three_geospiza_phased_neutral.vcf.gz

#bcftools view three_geospiza_phased_neutral.vcf -O z -o three_geospiza_phased_neutral.vcf.gz

#should have just done a for loop the way this worked out...
while read chr
do
    bcftools view -r $chr $GZVCF -O v -o ${chr}_three_geospiza_phased_neutral.vcf
    awk '{gsub(/^chr1A/,"199"); print}' ${chr}_three_geospiza_phased_neutral.vcf | awk '{gsub(/^chr4A/,"499"); print}' | awk '{gsub(/^chr/,""); print}'  > ${chr}_three_geospiza_phased_neutral_fixed.vcf #remove "chr" from the name
    $GEVA/geva_v1beta --vcf ${chr}_three_geospiza_phased_neutral_fixed.vcf --map $RECOMB/${chr}_recomb_rate.txt --out peak_${chr}
    awk 'NR%10000==1' peak_${chr}.marker.txt | grep -v "MarkerID"| cut -f 3 -d " " > ${chr}_trimmed_positions.txt #take a SNP every 10,000 lines
    $GEVA/geva_v1beta -i peak_${chr}.bin -o peak_${chr} --positions ${chr}_trimmed_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
done < $INTERVAL_LIST


# below written for peaks, will need to be adapted to select random vars
# while read chr start end locus
# do
#   ls ${locus}_positions.txt
#   ls $RECOMB/${chr}_recomb_rate.txt
#   $GEVA/geva_v1beta --vcf peak_${locus}_fixed.vcf --map $RECOMB/${chr}_recomb_rate.txt --out peak_${locus}
#   $GEVA/geva_v1beta -i peak_${locus}.bin -o peak_${locus} --positions ${locus}_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
# done < $PEAKS
#
#
# #need more for chr1A
# for locus in {2..10}
# do
#   awk '{gsub(/^chr1A/,"199"); print}' peak_${locus}.vcf > peak_${locus}_fixed.vcf
# done
#
# for locus in {2..10}
# do
#   chr=chr1A
#   $GEVA/geva_v1beta --vcf peak_${locus}_fixed.vcf --map $RECOMB/${chr}_recomb_rate.txt --out peak_${locus}
#   $GEVA/geva_v1beta -i peak_${locus}.bin -o peak_${locus} --positions ${locus}_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
# done < $PEAKS

#test runs I tried
# $GEVA/geva_v1beta -i chr1_test.bin -o chr1_test --positions chr1_test_positions.txt --Ne 11000 --mut 1.02e-8 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
# $GEVA/geva_v1beta -i chr1_test.bin -o chr1_test_ne40 --positions chr1_test_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 500 --maxDiscordant 500 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
# $GEVA/geva_v1beta -i chr1_test.bin -o chr1_test_ne40_max500 --positions chr1_test_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 500 --maxDiscordant 500 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
# $GEVA/geva_v1beta -i chr1_test.bin -o chr1_test_ne40_max1000 --positions chr1_test_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
