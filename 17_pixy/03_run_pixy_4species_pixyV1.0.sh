#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 8
#SBATCH -t 10-00:00:00
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

#NOTE SETUP TO RUN WITH ARRAY AND WHITE, YELLOW, OR OTHER
#FOR OTHER SHOULD ONLY RUN FOR PI
#SLURM_ARRAY_TASK_ID=20
POPS_FILE=four_geospiza_pixy.txt
VCF_DIR=/home/eenbody/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS/05_INCLUDE_INVARIANTS.bak/INV_PLUS_BIAL

#had to move this
OUTNAME=${POPS_FILE/_pixy.txt/}
mkdir -p $OUTNAME
cd $OUTNAME

#cut -f 1 $WORK_D/$POPS_FILE > ${OUTNAME}_samps.txt

#for SLURM_ARRAY_TASK_ID in {1..32}
#for SLURM_ARRAY_TASK_ID in 1
#do
  #VCF=$VCF_DIR/finches_sentieon_${SLURM_ARRAY_TASK_ID}_INV_BIALLELIC.vcf.gz
  #bcftools view -S ${OUTNAME}_samps.txt $VCF | bcftools annotate -x 'FORMAT,INFO' -O z -o finches_sentieon_${SLURM_ARRAY_TASK_ID}_INV_BIALLELIC_f4.vcf.gz
VCF=finches_sentieon_${SLURM_ARRAY_TASK_ID}_INV_BIALLELIC_f4.vcf.gz
tabix $VCF

#using slashes in zarr path and outfile prefix (i.e. must incl subdir) important
pixy --stats pi fst dxy \
  --vcf $VCF \
  --bypass_invariant_check yes \
  --n_cores 8 \
  --window_size 1000 \
  --populations $WORK_D/$POPS_FILE \
  --output_folder output_pixy1.0 \
  --output_prefix pixy_${OUTNAME}_${SLURM_ARRAY_TASK_ID}
#done
