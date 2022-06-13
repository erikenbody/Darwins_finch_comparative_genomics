#!/bin/bash
#SBATCH -A snic2018-8-65
#SBATCH -p core -n 7
#SBATCH -M snowy
#SBATCH -t 1-00:00:00
#SBATCH -J whatshap
#SBATCH -e whatshap_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o whatshap_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bioinfo-tools bcftools/1.10 samtools/1.10
#whatshap installed with pip
export PATH=$HOME/.local/bin:$PATH

WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/01_whatshap
cd $WORK_D

#31 chromosomes
INTERVAL_LIST=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/Cam_scaffolds.txt
INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
REF=/crex/proj/uppstore2017190/private/01_REFERENCE_DATA/Camarhynchus_parvulus_V1.0.fasta

#next step is to subset VCF for just parvulus
VCF=/crex/proj/uppstore2017190/private/04_VARIANTS/GATK_293/293.SNP.vcf.gz

####
SAMPLE=${1}

mkdir -p ${SNIC_TMP}/${SAMPLE}

bcftools view -r $INTERVAL -s $SAMPLE -O z -o ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.vcf.gz $VCF
bcftools index ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.vcf.gz

cd ${SNIC_TMP}/${SAMPLE}
scp $REF ${SNIC_TMP}/${SAMPLE}
scp $REF.fai ${SNIC_TMP}/${SAMPLE}

samtools view -b /crex/proj/uppstore2017190/private/03_ALIGNMENT/2017_Alignment/${SAMPLE}.MD.RG.bam $INTERVAL > ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.MD.RG.bam
samtools index ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.MD.RG.bam

whatshap phase -o ${SAMPLE}_${INTERVAL}.whatshap.vcf --reference=Camarhynchus_parvulus_V1.0.fasta ${SNIC_TMP}/${SAMPLE}/${SAMPLE}_${INTERVAL}.vcf.gz ${SAMPLE}_${INTERVAL}.MD.RG.bam

mkdir -p ${WORK_D}/${SAMPLE}
mv ${SAMPLE}_${INTERVAL}.whatshap.vcf ${WORK_D}/${SAMPLE}

bgzip ${WORK_D}/${SAMPLE}/${SAMPLE}_${INTERVAL}.whatshap.vcf
bcftools index ${WORK_D}/${SAMPLE}/${SAMPLE}_${INTERVAL}.whatshap.vcf.gz
