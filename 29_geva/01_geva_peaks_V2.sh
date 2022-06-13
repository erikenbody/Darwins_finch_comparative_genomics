#!/bin/bash
#SBATCH -A snic2021-5-251
#SBATCH -p core -n 8
#SBATCH -t 3-00:00:00
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

#bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.1 || ALT="*"' three_geospiza_phased.vcf.gz -O v -o three_geospiza_phased_var.vcf.gz

#geva input files
#use a bed file where 50kb extensions
PEAKS=$WORK_D/autosomes_for_ful_mag_species_lmm_PEAKS_extension.bed
GZVCF=$WORK_D/three_geospiza_phased_var.vcf.gz
RECOMB=$WORK_D/recombination_rate_per_chr_LDHELMET

mkdir -p peak_ages_V2
cd peak_ages_V2

# while read chr start end locus
# do
#     awk -v num=$locus 'NR == num' $PEAKS > tmp_peak.txt #franken code to ensure that tmp file has tab seperation
#     #cat tmp_peak.txt
#     #echo $locus
#     bedtools intersect -a $GZVCF -b tmp_peak.txt -wa -header > peak_${locus}.vcf
#     awk '{gsub(/^chr1A/,"199"); print}' peak_${locus}.vcf | awk '{gsub(/^chr4A/,"499"); print}' | awk '{gsub(/^chr/,""); print}' > peak_${locus}_fixed.vcf
#     rm tmp_peak.txt
# done < $PEAKS

while read chr start end locus
do
  ls ${locus}_positions.txt
  ls $RECOMB/${chr}_recomb_rate.txt

  grep -v "#" peak_${locus}_fixed.vcf | cut -f 2 > ${locus}_tmp.vcf #get list of positions in vcf
  grep -vf ${locus}_tmp.vcf $WORK_D/${locus}_positions.txt > ${locus}_missing_position.txt #find any positions not in vcf
  grep -vf  ${locus}_missing_position.txt $WORK_D/${locus}_positions.txt > ${locus}_positions_fixed.txt #remove the missing position

  sort ${locus}_positions_fixed.txt | head -n -1  > ${locus}_positions_fixed2.txt #get rid of last position because for these loci it crashes the run

  $GEVA/geva_v1beta --vcf peak_${locus}_fixed.vcf --map $RECOMB/${chr}_recomb_rate.txt --out peak_${locus}
  $GEVA/geva_v1beta -i peak_${locus}.bin -o peak_${locus} --positions ${locus}_positions_fixed2.txt --Ne 40000 --mut 1.02e-8 --maxConcordant 1000 --maxDiscordant 1000 --hmm $GEVA/hmm/hmm_initial_probs.txt $GEVA/hmm/hmm_emission_probs.txt
done < $PEAKS
