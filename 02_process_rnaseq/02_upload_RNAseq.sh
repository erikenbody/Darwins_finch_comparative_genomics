#!/bin/bash
#SBATCH -A snic2018-8-65
#SBATCH -M snowy
#SBATCH -p core -n 1
#SBATCH -t 4-00:00:00
#SBATCH -J ena_upload
#SBATCH -e ena_upload_%A_%a.err            # File to which STDERR will be written
#SBATCH -o ena_upload_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

MANIFEST=~/assem_finch_code/02_process_rnaseq/manifest_files/${SLURM_ARRAY_TASK_ID}_*.txt
WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/upload_RNAseq
FILE=$(basename $MANIFEST)
OUTPUT="${FILE%.*}"

mkdir $OUTPUT

java -Xmx4G -jar ${WORK_D}/webin-cli-1.8.4.jar -context reads -userName Webin-51919 -password "B,w'dbttO" -submit -inputDir $WORK_D -outputDir $OUTPUT -manifest $MANIFEST

#run like:
#sbatch --array=1-11 ~/assem_finch_code/02_process_rnaseq/02_upload_RNAseq.sh
