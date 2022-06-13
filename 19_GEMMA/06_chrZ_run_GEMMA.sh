#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 4
#SBATCH -M rackham
#SBATCH -t 00-06:00:00
#SBATCH -J gemma
#SBATCH -e gemma_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o gemma_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#Set up to run for PCs and also explory alternative covariates

#https://groups.google.com/forum/#!searchin/gemma-discussion/vcftools%7Csort:date/gemma-discussion/34FyII8oXUM/JAwdCu2mAQAJ

ml bioinfo-tools GEMMA/0.98.1 qctool/2.0.6 bcftools

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/GEMMA
GEMMA=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools
cd $WORK_D

#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
OUTNAME=for_ful_mag
CHR=chrZ
#############################################################


#PHENO=$WORK_D/for_ful_mag_pheno.txt
PHENO=$WORK_D/for_ful_mag_pheno_numerical.txt

##Comment again after chrZ and chr25 is run
VCF=finches_chrZ_324_shapeit4_phased_rh_GF.vcf.gz
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O z -o ${CHR}_fulldata_merge_ID.vcf.gz $VCF
# ##this makes a bit of a mess, but basically the problem is that qctools which makes the bimbim genotype chr_format
# ##outputs a snp ID as ID:chrom:pos. And the only ID thats easy to make is chrom_ID_pos_ref_alt
# ##could make it prettier by scripting my own qctool like thing, but why bother?
bcftools query -f '%ID\:%CHROM\:%POS, %POS\, %CHROM\n' ${CHR}_fulldata_merge_ID.vcf.gz > ${CHR}_fulldata_merge_map.txt
#
# ##MAKE SURE TO UNCOMMENT IF RE RUNNING THIS FROM SCRATCH
qctool -g ${CHR}_fulldata_merge_ID.vcf.gz -ofiletype bimbam_dosage -og ${OUTNAME}/${CHR}_${OUTNAME}_fulldata.geno

mkdir -p ${OUTNAME}/${CHR}_${OUTNAME}_gemma_out

$GEMMA/gemma -p $PHENO \
  -a ${CHR}_fulldata_merge_map.txt -g ${OUTNAME}/${CHR}_${OUTNAME}_fulldata.geno \
  -k $OUTNAME/autosomes_${OUTNAME}.noNA.relate.cXX.txt -lmm 2 -outdir ${OUTNAME}/${CHR}_${OUTNAME}_gemma_out \
  -n 2 \
  -o ${CHR}_${OUTNAME}_species_lmm_V2
