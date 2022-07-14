See the associated publication by Rubin and Enbody et al for the methods to construct this. Briefly, a LD-based 
linkage map was created for C. parvulus samples using LDhelmet. The resulting rho values were converted to centimorgan 
based on an Ne of 33,861. Methods are pasted below:

We computed recombination rate (rho) for 25 samples of C. parvulus from five islands using LDhelmet v1.10 (33). We 
additionally used three T. bicolor and five L. noctis samples from Barbados as outgroups. The phased haplotypes 
generated using SHAPEIT4 were converted to the required genotype format of LDhelmet using the vcftools –ldhelmet flag. 
In order to identify the ancestral state, we selected at each site the allele present in >60% of outgroup samples and 
following Singhal et al. (35) we assigned prior probabilities of 0.91 for the ancestral base and 0.03 of the remaining 
bases in order to account for the possibility of incorrect ancestral inference. At unresolved sites (e.g., missing 
data in outgroup samples), we used a prior of the stationary distribution of allele frequencies from the mutation rate 
matrix following the LDhelmet (33) method (also see Cam_parv_mutation_matrix.xlsx). Using genome-wide sites, we used 
the LDhelmet find_confs module with a window size of 50 SNPs to generate a haplotype configuration file. We next 
generated a likelihood lookup table using the table_gen command using the recommended rho grid values of (-r 0.0 0.1 
10.0 1.0 100.0) and a previously published
Watterson’s theta estimate of 0.00126(21). Last, we used the pade command to generate Padé coefficient tables using 
the same estimate of Watterson’s theta and the suggested 11 coefficients. We ran LDhelmet in 50 bp windows for 
10,000,000 iterations (discarding the first 200,000 as burn-in), and the following additional flags:
--max_lk_end 100 --prior_rate 0.05 -b 5.0
A previous study identified a block penalty of 5, as used here, as the ideal value for demonstrating fine scale 
recombination patterns (35). We extracted ρ/bp, representing the population scaled recombination rate, using 
post_to_text within LDhelmet. We summarized ρ by 20 kb sliding, non-overlapping windows (using the R package 
windowscanr, https://github.com/tavareshugo/WindowScanR), and report ρ/kb in comparison with diversity statistics. We 
additionally converted ρ to recombination fraction following r = ρ/2cNe (68) and Ne = 33,861 (21). Recombination 
fractions were converted to cM/Mb using the Morgan function r * 100 (for all plots in Fig. S2).
