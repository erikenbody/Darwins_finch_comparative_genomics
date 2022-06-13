#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 4
#SBATCH -t 4-00:00:00
#SBATCH -J gather
#SBATCH -e gather_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o gather_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.12

TOPDIR=/home/eenbody/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/07_BIALLELIC_ONLY/AC_FILTERED
cd $TOPDIR

# mkdir -p bgzipped

# for FILE in *.vcf.gz
# do
#   #bcftools view -O z -o bgzipped/${FILE} $FILE
#   #bcftools index bgzipped/${FILE}
#   ##did this then just reorganized into main directory
#   realpath ${FILE} >> list4merge.txt
# done

# bcftools concat --threads 4 -f list4merge.txt --naive -O z -o finches_sentieon_concat_SNPs_fil_setGT_PASS_BIALLELIC_ACf.vcf.gz
# bcftools index finches_sentieon_concat_SNPs_fil_setGT_PASS_BIALLELIC_ACf.vcf.gz

##make a list of invariant positions for external uses. It is somewhat practical to only include autosomes here. Chrz is interval 4
WORK_D=/home/eenbody/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/05_INCLUDE_INVARIANTS/INVARIANT_ONLY/
cd $WORK_D

mkdir -p pos_files

for FILE in *.vcf.gz
do
  SAMPBASE=$(basename $FILE)
  SAMPLE=${SAMPBASE/.vcf.gz}
  bcftools view -s 00Esp1 -O z -o pos_files/${SAMPLE}_00Esp1.vcf.gz $FILE
  realpath pos_files/${SAMPLE}_00Esp1.vcf.gz | grep -v "sentieon_4_INVARIANT" >> list4merge2.txt
done
bcftools concat --threads 4 -f list4merge2.txt --naive -O z -o finches_sentieon_autosomes_INVARIANT.vcf.gz
