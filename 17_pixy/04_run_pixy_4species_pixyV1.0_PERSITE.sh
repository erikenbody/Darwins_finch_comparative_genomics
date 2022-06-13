#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 8
#SBATCH -t 5-00:00:00
#SBATCH -J pixy
#SBATCH -e pixy_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o pixy_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml conda/latest bioinfo-tools bcftools/1.10
source conda_init.sh

conda activate pixy

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/pixy
cd $WORK_D

POPS_FILE=four_geospiza_pixy.txt
VCF_DIR=/home/eenbody/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/05_INCLUDE_INVARIANTS.bak/INV_PLUS_BIAL

#had to move this
OUTNAME=${POPS_FILE/_pixy.txt/}
mkdir -p $OUTNAME
cd $OUTNAME

VCF=finches_sentieon_${SLURM_ARRAY_TASK_ID}_INV_BIALLELIC_f4.vcf.gz

#ONLY FST

pixy --stats fst \
  --vcf $VCF \
  --bypass_invariant_check yes \
  --n_cores 8 \
  --window_size 1 \
  --populations $WORK_D/$POPS_FILE \
  --output_folder output_pixy1.0_PERSITE \
  --output_prefix pixy_PERSITE_${OUTNAME}_${SLURM_ARRAY_TASK_ID}
