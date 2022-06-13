#!/bin/bash
#SBATCH -A snic2018-8-65
#SBATCH -p core -n 4
#SBATCH -M rackham
#SBATCH -t 2-00:00:00
#SBATCH -J whatshap
#SBATCH -e whatshap_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o whatshap_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10

WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/01_ancestral_state
cd $WORK_D
SI=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/tools/shapeit4/bin
INTERVAL_LIST=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/Cam_scaffolds.txt
INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

module load ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10
$SI/shapeit4 --input barbados.whatshap.vcf.gz \
          --use-PS 0.0001 \
          --output barbados.shapeit4.phased_${INTERVAL}.bcf \
          --thread 4 \
          --log phased.log \
          --region $INTERVAL

#sbatch --array=1-31 ~/afc/10_LDhelmet/01.2_anc_state_shapeit.sh
