#!/bin/bash
#SBATCH -A snic2020-5-685
#SBATCH -p core -n 4
#SBATCH -M rackham #snowy was slow for wait time
#SBATCH -t 1-0:00:00
#SBATCH -J fasttree
#SBATCH -e fasttree_%J_%A_%a.err            # File to which STDERR will be written
#SBATCH -o fasttree_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=erik.enbody@gmail.com

ml bcftools/1.10 BEDTools/2.29.2
ml python/2.7.15 FastTree/2.1.8

WORK_D=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/phylogeny
cd $WORK_D
DROP_BED=autosomes_for_ful_mag_species_lmm_PEAKS.bed
VCF2PHYLIP=/crex/proj/snic2020-2-19/private/wagtail/users/erikenbody/X_larks/02_phylogenetics/vcf2phylip/

GZVCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/GEMMA/finches_324_shapeit4_phased_rh_GF_prune.vcf.gz
bedtools intersect -v -a $GZVCF -b $DROP_BED -wa -header > finches_324_shapeit4_phased_rh_GF_prune_autosomes_neutral.vcf
python $VCF2PHYLIP/vcf2phylip.py --resolve-IUPAC --input finches_324_shapeit4_phased_rh_GF_prune_autosomes_neutral.vcf --fasta
FastTree -nt -log finch_fasttree_log_file -gtr -fastest finches_324_shapeit4_phased_rh_GF_prune_autosomes_neutral.min4.fasta > finches_324_shapeit4_phased_rh_GF_prune_autosomes_neutral.min4.tre

##no pruning.
bedtools intersect -v -a $GZVCF -b $DROP_BED -wa -header > finches_324_shapeit4_phased_rh_GFautosomes_neutral.vcf
python $VCF2PHYLIP/vcf2phylip.py --resolve-IUPAC --input finches_324_shapeit4_phased_rh_GFautosomes_neutral.vcf --fasta
FastTree -nt -log finch_fasttree_log_file -gtr -fastest finches_324_shapeit4_phased_rh_GFautosomes_neutral.min4.fasta > finches_324_shapeit4_phased_rh_GFautosomes_neutral.min4.tre
#
#
# ####do the same but only for ROI
GZVCF=/crex/proj/snic2020-2-19/private/darwins_finches/users/erikenbody/Finch_asm/GEMMA/finches_324_shapeit4_phased_rh_GF.vcf.gz
#
bedtools intersect -a $GZVCF -b $DROP_BED -wa -header > finches_324_shapeit4_phased_rh_GFautosomes_peaks.vcf
python $VCF2PHYLIP/vcf2phylip.py --resolve-IUPAC --input finches_324_shapeit4_phased_rh_GFautosomes_peaks.vcf --fasta
FastTree -nt -log finch_fasttree_log_file -gtr -fastest finches_324_shapeit4_phased_rh_GFautosomes_peaks.min4.fasta > finches_324_shapeit4_phased_rh_GFautosomes_peaks.min4.tre
