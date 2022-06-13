#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 5
#SBATCH -t 4-00:00:00
#SBATCH -J gemma-wrapper
#SBATCH -e gemma-wrapper_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o gemma-wrapper_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#https://groups.google.com/forum/#!searchin/gemma-discussion/vcftools%7Csort:date/gemma-discussion/34FyII8oXUM/JAwdCu2mAQAJ

PATH=$PATH:/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools

ml bioinfo-tools GEMMA/0.98.1 qctool/2.0.6 bcftools/1.10

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/GEMMA
GEMMA=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_lowpass/tools
cd $WORK_D

#####JUST CHANGE THIS TO RE RUN WITH DIFFERENT SAMPLES#######
OUTNAME=for_ful_mag
CHR=autosomes
#############################################################

PHENO=$WORK_D/for_ful_mag_pheno_numerical_TAB_SEP.txt

mkdir -p ${OUTNAME}/${CHR}_${OUTNAME}_gemma_out

GENO=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/GEMMA/for_ful_mag/autosomes_for_ful_mag_noNA.geno

# gemma-wrapper --cache-dir ./${OUTNAME}/${OUTNAME}_gemma_cache --json -- \
#     -g $GENO \
#     -p $PHENO \
#     -gk \
#     -debug > ${OUTNAME}/${OUTNAME}_K.json

# echo "k done"

NUM_PHENO=2
cut -f $NUM_PHENO $PHENO > ${OUTNAME}/${OUTNAME}_${NUM_PHENO}.txt

gemma-wrapper --input ${OUTNAME}/${OUTNAME}_pruned_K.json --permutate 100 --permute-phenotype ${OUTNAME}/${OUTNAME}_${NUM_PHENO}.txt -- \
   -g $OUTNAME/${CHR}_${OUTNAME}.geno \
   -a autosomes_fulldata_merge_map.txt \
   -lmm 2 \
   -debug > ${OUTNAME}/${OUTNAME}_${NUM_PHENO}_GWA_100_lmm2.txt
