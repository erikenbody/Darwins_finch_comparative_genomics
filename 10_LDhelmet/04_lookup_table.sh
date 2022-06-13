#!/bin/bash
#SBATCH -A snic2018-8-63
#SBATCH -p node -n 20
#SBATCH -M rackham
#SBATCH -t 10-00:00:00
#SBATCH -J lookup_tables
#SBATCH -e lookup_tables_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o lookup_tables_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

#this only took a day or two to run, mostly by the pade calculations

ml LDhelmet/1.10
WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/04_lookup_tables
cd $WORK_D
SAMPLE_DIR=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/LDhelmet/03_fasta_format

realpath $SAMPLE_DIR/*.snps | grep -v "chrZ" > ldhelmet_input.txt

#theta set based on 2015 finch paper nucleotide diversity

ldhelmet find_confs --num_threads 20 -w 50 -o parv.conf $(<ldhelmet_input.txt)
echo "find_confs done"
ldhelmet table_gen --num_threads 20 -c parv.conf -t 0.00126 -r 0.0 0.1 10.0 1.0 100.0 -o parv.lk
echo "lk tabls done"

ldhelmet pade --num_threads 20 -c parv.conf -t 0.00126 -x 11 -o parv.pade
echo "pade done"
#sbatch -p devel -t 00:10:00 ~/afc/10_LDhelmet/04_lookup_table.sh
