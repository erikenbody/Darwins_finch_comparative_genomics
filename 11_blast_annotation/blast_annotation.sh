#!/bin/bash -l
#SBATCH -A snic2018-8-65
#SBATCH -M snowy
#SBATCH -p node -n 6
#SBATCH -t 2-0:00:00
#SBATCH -J blastp
#SBATCH -e blastp_%A_%a.err            # File to which STDERR will be written
#SBATCH -o blastp_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

module load bioinfo-tools maker/3.01.2-beta blast/2.7.1+

WORK_D=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/blast_annotation
cd $WORK_D

makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot_ONLY -dbtype prot -title uniprot_sprot_ONLY

#reformat names for Annie later. Specifically need to use transcript names
sed 's/^.*\(transcript:\)/\1/' Camarhynchus_parvulus.Camarhynchus_parvulus_V1.1.pep.all.fa | sed 's/transcript:/>transcript:/g' | sed 's/\.1 gene.*//g' > Camarhynchus_parvulus.Camarhynchus_parvulus_V1.1_transcript_names.pep.all.fa

blastp -db $WORK_D/uniprot_sprot_ONLY -query Camarhynchus_parvulus.Camarhynchus_parvulus_V1.1_transcript_names.pep.all.fa -out Camarhynchus_parvulus_V1.1_blast_swissprot.out -evalue .000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads 6

ANNIE=/crex/proj/uppstore2019015/private/DarwinsFinches/finch_assembly/tools/annie/genomeannotation-annie-4bb3980/

#not needed because of edits above
#sed 's/\.1//g' Camarhynchus_parvulus_V1.1_blast_swissprot.out > Camarhynchus_parvulus_V1.1_blast_swissprot_4Annie.out

python3 $ANNIE/annie.py -b Camarhynchus_parvulus_V1.1_blast_swissprot.out \
        -g Camarhynchus_parvulus.Camarhynchus_parvulus_V1.1.100.gff3 \
        -db uniprot_sprot.fasta \
        --output Camarhynchus_parvulus_V1.1.annie
