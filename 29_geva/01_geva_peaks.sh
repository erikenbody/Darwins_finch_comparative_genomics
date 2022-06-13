#!/bin/bash
#SBATCH -A snic2021-5-251
#SBATCH -p core -n 8
#SBATCH -t 1-00:00:00
#SBATCH -J geva
#SBATCH -e geva_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o geva_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.14 BEDTools/2.29.2

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/geva
cd $WORK_D

GEVA=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/tools/geva

#run once to subset phased vcf file
#SAMPS=3_geospiza_sp.txt
#VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/phasing/finches_324_shapeit4_phased.bcf
#bcftools view -S $SAMPS -O z -o three_geospiza_phased.vcf.gz $VCF

#geva input files
PEAKS=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/phylogeny/autosomes_for_ful_mag_species_lmm_PEAKS.bed
GZVCF=three_geospiza_phased.vcf.gz
RECOMB=recombination_rate_per_chr_LDHELMET

#should have just done a for loop the way this worked out...
# while read chr start end locus
# do
#     awk -v num=$locus 'NR == num' $PEAKS > tmp_peak.txt #franken code to ensure that tmp file has tab seperation
#     #cat tmp_peak.txt
#     #echo $locus
#     bedtools intersect -a $GZVCF -b tmp_peak.txt -wa -header > peak_${locus}.vcf
#     awk '{gsub(/^chr/,""); print}' peak_${locus}.vcf > peak_${locus}_fixed.vcf
#     rm tmp_peak.txt
# done < $PEAKS
#
#
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

for locus in 10 20 21 28 9
do
  grep -v "#" peak_${locus}_fixed.vcf | cut -f 2 > ${locus}_tmp.vcf #get list of positions in vcf
  grep -vf ${locus}_tmp.vcf ${locus}_positions.txt > ${locus}_missing_position.txt #find any positions not in vcf
  grep -vf  ${locus}_missing_position.txt ${locus}_positions.txt > ${locus}_positions_fixed.txt #remove the missing position

  sort ${locus}_positions_fixed.txt | head -n -1  > ${locus}_positions_fixed2.txt #get rid of last position because for these loci it crashes the run

  $GEVA/geva_v1beta -t 8 -i peak_${locus}.bin -o peak_${locus} --positions ${locus}_positions_fixed2.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
done

#how i found that the last position causes an erro
# while read POS
# do
#   echo $POS
#   $GEVA/geva_v1beta -t 8 -i peak_${locus}.bin -o ${POS}_${locus} --position $POS --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
# done < ${locus}_positions_fixed.txt
#  # problematic _positions_fixed
# $GEVA/geva_v1beta -t 8 -i peak_${locus}.bin -o 50899028_${locus} --position 50899028 --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt


#test runs I tried
# $GEVA/geva_v1beta -i chr1_test.bin -o chr1_test --positions chr1_test_positions.txt --Ne 11000 --mut 1.02e-8 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
# $GEVA/geva_v1beta -i chr1_test.bin -o chr1_test_ne40 --positions chr1_test_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 500 --maxDiscordant 500 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
# $GEVA/geva_v1beta -i chr1_test.bin -o chr1_test_ne40_max500 --positions chr1_test_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 500 --maxDiscordant 500 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
# $GEVA/geva_v1beta -i chr1_test.bin -o chr1_test_ne40_max1000 --positions chr1_test_positions.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
