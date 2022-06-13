#!/bin/bash
#SBATCH -A snic2018-8-65
#SBATCH -p core -n 1 #was full node mem 256gb for whatshap
#SBATCH -M rackham
#SBATCH -t 10-00:00:00
#SBATCH -J whatshap
#SBATCH -e whatshap_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o whatshap_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10

#this script would run much much faster if it were run per chromosome

WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/01_ancestral_state
cd $WORK_D
SI=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/tools/shapeit4/bin

REF=/crex/proj/uppstore2017190/private/01_REFERENCE_DATA/Camarhynchus_parvulus_V1.0.fasta

#next step is to subset VCF for just outgroup taxa
VCF=/crex/proj/uppstore2017190/private/04_VARIANTS/GATK_293/293.SNP.vcf.gz
bcftools query -l $VCF | grep "Barb" > barb_names.txt
realpath /crex/proj/uppstore2017190/private/03_ALIGNMENT/2017_Alignment/*Barb*bam > barb_bamlist.txt

bcftools view -S barb_names.txt -o barbados.vcf $VCF
bcftools index $VCF

#installed whatshap with pip
export PATH=$HOME/.local/bin:$PATH
whatshap phase -o barbados.whatshap.vcf --reference=$REF barbados.vcf $(<barb_bamlist.txt)

bgzip barbados.whatshap.vcf
bcftools index barbados.whatshap.vcf.gz
