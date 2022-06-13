#!/bin/bash
#SBATCH -A snic2021-22-187
#SBATCH -p core -n 2
#SBATCH -t 1-00:00:00
#SBATCH -J tajd
#SBATCH -e tajd_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o tajd_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bioinfo-tools vcftools

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/tajimas_d
cd $WORK_D

VCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/GEMMA/finches_324_shapeit4_phased_rh_GF.vcf.gz

# vcftools --gzvcf $VCF --keep fortis_samples.txt --TajimaD 20000 --out fortis_tajd
# vcftools --gzvcf $VCF --keep fuliginosa_samples.txt --TajimaD 20000 --out fuliginosa_tajd
# vcftools --gzvcf $VCF --keep magnirostris_samples.txt --TajimaD 20000 --out magnirostris_tajd

#with invariant

mkdir -p vcftools_with_invariants
vcftools --gzvcf vcfkit/fortis_${SLURM_ARRAY_TASK_ID}.vcf.gz --TajimaD 20000 --out fortis_V2_tajd_${SLURM_ARRAY_TASK_ID}
vcftools --gzvcf vcfkit/fuliginosa_${SLURM_ARRAY_TASK_ID}.vcf.gz --TajimaD 20000 --out fuliginosa_V2_tajd_${SLURM_ARRAY_TASK_ID}
vcftools --gzvcf vcfkit/magnirostris_${SLURM_ARRAY_TASK_ID}.vcf.gz --TajimaD 20000 --out magnirostris_V2_tajd_${SLURM_ARRAY_TASK_ID}
