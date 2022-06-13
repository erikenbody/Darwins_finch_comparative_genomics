#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 4 #was 4
#SBATCH -t 2-00:00:00
#SBATCH -J set_filters
#SBATCH -e set_filters_%J_%A_%a.err
#SBATCH -o set_filters_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bioinfo-tools samtools vcftools bcftools
module load GATK/4.1.4.1

TOPDIR=/home/eenbody/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS
cd $TOPDIR

REF=/crex/proj/uppstore2017190/private/01_REFERENCE_DATA/Camarhynchus_parvulus_V1.0.fasta
INTERVAL_FILE=/crex/proj/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/Cam_scaffolds_corrected.txt
INTERVAL=`cat $INTERVAL_FILE | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
#1-32, 33 is unplaced

VCF=finches_sentieon_415_rh.vcf.gz

##updating to run with V2 and all one script. Dont need 03 select now
gatk --java-options "-Xmx29g" SelectVariants  \
    -R $REF \
    -V $VCF \
    --select-type-to-include SNP \
    --select-type-to-include NO_VARIATION \
    -O finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV.vcf.gz \
    --intervals $INTERVAL

#set filter expressions

gatk --java-options "-Xmx29g" VariantFiltration  \
    -R $REF \
    -V finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV.vcf.gz \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0))" \
    --filter-name "RPRS8" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('QD') && QD < 2.0))" \
    --filter-name "QD2" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('FS') && FS > 60.0))" \
    --filter-name "FS60" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('SOR') && SOR > 3.0))" \
    --filter-name "SOR3" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('MQ') && MQ < 40.0))" \
    --filter-name "MQ40" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
    --filter-name "MQ12.5" \
    -G-filter "vc.isSNP() && DP < 2" \
    -G-filter-name "gtDP1" \
    -G-filter "vc.isSNP() && DP > 100" \
    -G-filter-name "gtDP150" \
    -G-filter "vc.isSNP() && GQ < 20" \
    -G-filter-name "gtGQ10" \
    -O finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil.vcf.gz

#set filtered GT to no call
gatk --java-options "-Xmx29g" SelectVariants  \
    -R $REF \
    -V finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil.vcf.gz \
    --set-filtered-gt-to-nocall \
    -O finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT.vcf.gz

mkdir -p 03_FILTERS_SET
vcftools --gzvcf finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT.vcf.gz --FILTER-summary --out  03_FILTERS_SET/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT.vcf.summary

mkdir -p 05_INCLUDE_INVARIANTS
#filtered with invariants. new line here
#bcftools view --threads 4 -f PASS -e 'ALT="*" | TYPE~"indel" | ref="N"' finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT.vcf.gz | bcftools annotate -x '^FORMAT/AN,FORMAT/AC,^INFO/AD,INFO/DP' -O z -o 05_INCLUDE_INVARIANTS/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT_PASS.vcf.gz
bcftools view --threads 4 -f PASS -e 'ALT="*" | TYPE~"indel" | ref="N"' finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT.vcf.gz | bcftools annotate -x '^INFO/AN,INFO/AC,^FORMAT/GT,FORMAT/AD,FORMAT/DP' -O z -o 05_INCLUDE_INVARIANTS/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT_PASS.vcf.gz
tabix 05_INCLUDE_INVARIANTS/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT_PASS.vcf.gz

mkdir -p 06_SNPS_ONLY
#filtered SNPs only
bcftools view --threads 4 -v snps -O z -o 06_SNPS_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_fil_setGT_PASS.vcf.gz 05_INCLUDE_INVARIANTS/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT_PASS.vcf.gz
tabix 06_SNPS_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_fil_setGT_PASS.vcf.gz

mkdir -p 07_BIALLELIC_ONLY
#filtered biallelic SNPs only
bcftools view --threads 4 -v snps -m2 -M2 -O z -o 07_BIALLELIC_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_fil_setGT_PASS_BIALLELIC.vcf.gz 06_SNPS_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_fil_setGT_PASS.vcf.gz
tabix 07_BIALLELIC_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_fil_setGT_PASS_BIALLELIC.vcf.gz

#apply a filter to remove invariant lingering and also singletons
mkdir -p 07_BIALLELIC_ONLY/AC_FILTERED
bcftools view --threads 4 -e 'AC<=1 || AC==AN' -O z -o 07_BIALLELIC_ONLY/AC_FILTERED/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_fil_setGT_PASS_BIALLELIC_ACf.vcf.gz 07_BIALLELIC_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_fil_setGT_PASS_BIALLELIC.vcf.gz

VCF=05_INCLUDE_INVARIANTS/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_INV_fil_setGT_PASS.vcf.gz
mkdir -p 05_INCLUDE_INVARIANTS/INVARIANT_ONLY

#require that 50% of samples were non missing to include invariant site
#note that later I found that this does not exclude triallelic sites
vcftools --gzvcf $VCF \
--max-maf 0 \
--max-missing 0.50 \
--max-alleles 1 \
--recode --stdout | bgzip -c > 05_INCLUDE_INVARIANTS/INVARIANT_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_INVARIANT.vcf.gz

tabix 05_INCLUDE_INVARIANTS/INVARIANT_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_INVARIANT.vcf.gz

VCF=05_INCLUDE_INVARIANTS/INVARIANT_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_INVARIANT.vcf.gz
VCF_SNPS=07_BIALLELIC_ONLY/finches_sentieon_${SLURM_ARRAY_TASK_ID}_SNPs_fil_setGT_PASS_BIALLELIC.vcf.gz

mkdir -p 05_INCLUDE_INVARIANTS/INV_PLUS_BIAL

bcftools concat \
--allow-overlaps \
--rm-dups all \
$VCF_SNPS $VCF \
-O z -o 05_INCLUDE_INVARIANTS/INV_PLUS_BIAL/finches_sentieon_${SLURM_ARRAY_TASK_ID}_INV_BIALLELIC.vcf.gz

tabix 05_INCLUDE_INVARIANTS/INV_PLUS_BIAL/finches_sentieon_${SLURM_ARRAY_TASK_ID}_INV_BIALLELIC.vcf.gz
