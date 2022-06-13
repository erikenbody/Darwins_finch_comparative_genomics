#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -p core -n 4
#SBATCH -M rackham
#SBATCH -t 0-05:00:00
#SBATCH -J shapeit4
#SBATCH -e shapeit4_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o shapeit4_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10

WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/02_shapeit4
cd $WORK_D

SI=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/tools/shapeit4/bin
INTERVAL_LIST=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/Cam_scaffolds.txt
INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
SAMPLE_DIR=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/01_whatshap

ls $SAMPLE_DIR/parvulus_*/*${INTERVAL}.whatshap.vcf.gz > ${INTERVAL}.list

bcftools merge -l ${INTERVAL}.list -O z -o ${INTERVAL}.whatshap.vcf.gz
bcftools index ${INTERVAL}.whatshap.vcf.gz

module load ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10
$SI/shapeit4 --input ${INTERVAL}.whatshap.vcf.gz \
          --use-PS 0.0001 \
          --output ${INTERVAL}.parvulus.phased.vcf.gz \
          --thread 4 \
          --region $INTERVAL \
          --log ${INTERVAL}.phasing.log

#sbatch --array=2-29 ~/afc/10_LDhelmet/02_parv_shapeit.sh
#sbatch --array=1 ~/afc/10_LDhelmet/02_parv_shapeit.sh
