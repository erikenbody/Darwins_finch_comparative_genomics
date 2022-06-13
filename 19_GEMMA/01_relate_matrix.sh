#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 1 #dont need much to run this on the pruned dataset. interactive is fine
#SBATCH -t 0-05:00:00
#SBATCH -J related
#SBATCH -e related_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o related_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#https://groups.google.com/forum/#!searchin/gemma-discussion/vcftools%7Csort:date/gemma-discussion/34FyII8oXUM/JAwdCu2mAQAJ

ml bioinfo-tools GEMMA/0.98.1 qctool/2.0.6

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/GEMMA
GEMMA=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools
cd $WORK_D

BCFTOOLS=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/bcftools_devel/bcftools
export BCFTOOLS_PLUGINS=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools/bcftools_devel/bcftools/plugins

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/phasing/finches_324_shapeit4_phased_rh.bcf

#note this drops chrZ. Then removes any sites no longer variable
$BCFTOOLS/bcftools view -t ^chrZ -S for_ful_mag_samples.txt $VCF | $BCFTOOLS/bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.1 || ALT="*"' -O v -o finches_324_shapeit4_phased_rh_GF.vcf.gz
$BCFTOOLS/bcftools +prune -n 1 -w 20000bp --nsites-per-win-mode rand --random-seed 42 -O v -o finches_324_shapeit4_phased_rh_GF_prune.vcf.gz finches_324_shapeit4_phased_rh_GF.vcf.gz
$BCFTOOLS/bcftools index finches_324_shapeit4_phased_rh_GF_prune.vcf.gz

VCF=finches_324_shapeit4_phased_rh_GF.vcf.gz

CHR=autosomes
#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
OUTNAME=for_ful_mag
#############################################################

PHENO=$WORK_D/for_ful_mag_pheno.txt

mkdir -p $OUTNAME

##formatting of map and SNP id only ever needs to be run once
$BCFTOOLS/bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O z -o ${CHR}_merge_ID.vcf.gz $VCF
##this makes a bit of a mess, but basically the problem is that qctools which makes the bimbim genotype chr_format
##outputs a snp ID as ID:chrom:pos. And the only ID thats easy to make is chrom_ID_pos_ref_alt
##could make it prettier by scripting my own qctool like thing, but why bother?
$BCFTOOLS/bcftools query -f '%ID\:%CHROM\:%POS, %POS\, %CHROM\n' ${CHR}_merge_ID.vcf.gz > ${CHR}_merge_map.txt

##Can start here if Daphne birds

qctool -g ${CHR}_merge_ID.vcf.gz -ofiletype bimbam_dosage -og $OUTNAME/${CHR}_${OUTNAME}.geno

##create relatedness matrix (centered is best reccomendation for starting off in the documentation)

grep -v "NA" $OUTNAME/${CHR}_${OUTNAME}.geno  > $OUTNAME/${CHR}_${OUTNAME}_noNA.geno

$GEMMA/gemma -g $OUTNAME/${CHR}_${OUTNAME}_noNA.geno -p $PHENO -gk 1 -o ${CHR}_${OUTNAME}.noNA.relate
mv output/* $OUTNAME
