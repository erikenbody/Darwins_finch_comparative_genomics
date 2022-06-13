#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 10
#SBATCH -M rackham
#SBATCH -t 5-00:00:00
#SBATCH -J sentieon
#SBATCH -e sentieon_%A_%a.err            # File to which STDERR will be written
#SBATCH -o sentieon_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml python/2.7.15

TOPDIR=/home/eenbody/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS
cd $TOPDIR

REFERENCE=/crex/proj/uppstore2017190/private/01_REFERENCE_DATA/Camarhynchus_parvulus_V1.0.fasta

FILENAME=`cat TH-2647_Feb21_new_fastq.txt |  awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
SAMPBASE=$(basename $FILENAME)
SAMPLE=${SAMPBASE%_S*}
SAMPLE=${SAMPLE/TH-2647-/}

SAMPLE_DIR=$(dirname $FILENAME)

R1=$SAMPLE_DIR/*${SAMPLE}*R1*fastq.gz
R2=$SAMPLE_DIR/*${SAMPLE}*R2*fastq.gz
echo $SAMPLE

if [ -f ${SAMPLE}/${SAMPLE}.sort.bam ]; then echo "bam exists"; else sh ~/afc/16_sentieon_invariants/sentieon_finches_2021.sh $R1 $R2 $REFERENCE $SAMPLE $SAMPBASE; fi

##quickly make sure these files truely exist
while read FILENAME
do
  SAMPBASE=$(basename $FILENAME)
  SAMPLE=${SAMPBASE%_S*}
  SAMPLE=${SAMPLE/TH-2647-/}
  ll $SAMPLE/${SAMPLE}.g.vcf.gz >> new_gvcfs.txt
done < TH-2647_Feb21_new_fastq.txt
