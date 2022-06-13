#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p node
#SBATCH -t 7-00:00:00
#SBATCH -J genotyper
#SBATCH -e genotyper_%A_%a.err            # File to which STDERR will be written
#SBATCH -o genotyper_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml python/2.7.15

TOPDIR=/home/eenbody/snic2020-2-19/private/darwins_finches/variants/INCLUDE_INVARIANTS
cd $TOPDIR

REFERENCE=/crex/proj/uppstore2017190/private/01_REFERENCE_DATA/Camarhynchus_parvulus_V1.0.fasta

#sentieon options
#export SENTIEON_LICENSE=/crex/proj/uppstore2019097/nobackup/enbody_wagtails_working/tools/sentieon/Uppsala_University_eval.lic
export SENTIEON_LICENSE=/domus/h1/eenbody/Uppsala_cluster3.lic

SENTIEON_INSTALL_DIR=/crex/proj/uppstore2019097/nobackup/enbody_wagtails_working/tools/sentieon/sentieon-genomics-201911
SENTIEON_TMPDIR=$SNIC_TMP

#make list of variants. But only do once
# for GVCF in */*.g.vcf.gz
# do
#   LST=gvcf_count_V3.txt
#   ls $GVCF >> gvcf_count_V3.txt
#   echo -ne " -v $GVCF" >> gvcf_list_V3.txt
# done
#
# wc -l gvcf_count_V3.txt
# wc -l gvcf_count_V2.txt

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REFERENCE -t 20 --algo GVCFtyper --emit_mode CONFIDENT finches_sentieon_416.vcf.gz $(<gvcf_list_V3.txt)



#$SENTIEON_INSTALL_DIR/bin/sentieon licsrvr --start --log Mar3_sent.log Uppsala_cluster3.lic
