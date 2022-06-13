#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -p core -n 5
#SBATCH -M rackham
#SBATCH -t 10-00:00:00
#SBATCH -J ldhelmet
#SBATCH -e ldhelmet_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o ldhelmet_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bioinfo-tools LDhelmet/1.10

#modules to load to run as custom install, which didnt improve multi threading like I hoped so not worth it
#module load bioinfo-tools ABINIT/8.10.3 GCCcore/8.3.0 bcftools/1.10
#export LD_LIBRARY_PATH=/crex/proj/uppstore2017190/b2012111_nobackup/private/Erik/finches/tools/boost_1_73_0/lib:$LD_LIBRARY_PATH
#LDHELMET_DIR=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/tools/LDhelmet_v1.10

WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/05_run_LDHelmet
cd $WORK_D

#intervals for running per chromosome
INTERVAL_LIST=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/Cam_scaffolds.txt
INTERVAL=`cat $INTERVAL_LIST | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
SAMPLE_DIR=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/03_fasta_format

#LDhelmet paramenter input
MUTATION_MATRIX=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/01_ancestral_state/normalized_parvulus_mutation.matrix
ANC_DIR=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/01_ancestral_state
PADE=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/04_lookup_tables/parv.pade
LK=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/04_lookup_tables/parv.lk

#the ancestral state must be determined for every SNP. For some, it cannot be determined. In the Przeworski paper, they fill it in:
#"by providing prior probabilities for each nucleotide that were equal to their stationary frequency, as estimated from empirical frequencies in the genome and the mutation matrix "
#So, I used the values from the diagonal in the mutation matrix for unknown positions

#mats suggestion
#maf < 0.1
#200,000 burn in
#-n 10000000

#make sure to update this with appropriate values
#cut -f 3-6 $ANC_DIR/${INTERVAL}.anc | sed 's/0\t0\t0\t0/0.358\t0.002\t0.000\t0.360/g' > ${INTERVAL}_TMP.anc
#average of columns in mutation matrix
cut -f 3-6 $ANC_DIR/${INTERVAL}.anc | sed 's/0\t0\t0\t0/0.320\t0.180\t0.180\t0.321/g' > ${INTERVAL}_TMP.anc
cut -f 2 $ANC_DIR/${INTERVAL}.anc > ${INTERVAL}_TMP.pos
paste ${INTERVAL}_TMP.pos ${INTERVAL}_TMP.anc > ${INTERVAL}_recalibrate.anc

#ldhelmet rjmcmc --num_threads 5 -w 50 -l $LK -p $PADE -b 5.0 --snps_file $SAMPLE_DIR/${INTERVAL}.parvulus.ldhelmet.snps --pos_file $SAMPLE_DIR/${INTERVAL}.parvulus.ldhelmet.pos -m $MUTATION_MATRIX -a ${INTERVAL}_recalibrate.anc --burn_in 100000 -n 1000000 -o ${INTERVAL}_parvulus.post

#re running more burn in and nor n, also adding max_lk_end and prior_rate to match singhal science 2015 paper
ldhelmet rjmcmc --num_threads 5 -w 50 -l $LK -p $PADE -b 5.0 --snps_file $SAMPLE_DIR/${INTERVAL}.parvulus.ldhelmet.snps --pos_file $SAMPLE_DIR/${INTERVAL}.parvulus.ldhelmet.pos -m $MUTATION_MATRIX -a ${INTERVAL}_recalibrate.anc --burn_in 200000 -n 10000000 -o ${INTERVAL}_parvulus.post --max_lk_end 100 --prior_rate 0.05

ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.0975 -o ${INTERVAL}_parvulus_ldhelmet.txt ${INTERVAL}_parvulus.post

#sbatch --array=1-3,5-31 --dependency=afterok:14756844 ~/afc/10_LDhelmet/05_run_LDHelmet.sh
