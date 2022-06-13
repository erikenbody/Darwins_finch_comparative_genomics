# Load the necessary libraries
library(ape) # For plotting phylogenetic trees
library(grid) # Used to plot lines between plot panels
library(ggtree)



#chZ<-read.tree("output/svdquartets/sweep/allquartets/wagtails_214_Pseudosweep.SNPs.fil.PASS.SPANRMV.GENOFILT.sub.thin.NEWICK.tre")
#auto<-read.tree("output/svdquartets/auto/wagtails_214_Pseudochr8.SNPs.fil.PASS.SPANRMV.GENOFILT.sub.thin.OUT.tre")

#auto<-read.tree("data/newick_trees/rooted_finches_324_shapeit4_phased_rh_maf0.15.recode_no_region_selection.for.poptree.txt.for.poptree.d3.b1.tre")
#sweep<-read.tree("data/newick_trees/rooted_Top100SNPs_gwas_loci_d3_100.nj.tre")

auto1 <- read.tree("data/newick_trees/science_hybrid_paper_species_phylogeny_TRIMMED_TIME_RENAMED.tre")

auto1$tip.label <- gsub("Platyspiza_","P.",gsub("Tiaris_","T.",gsub("Certhidea_","C.",gsub("Loxigilla_","L.",gsub("Pinaroloxias_","P.",gsub("Camarhynchus_","C.",gsub("Geospiza_","G.",auto1$tip.label)))))))
write.tree(auto1, "data/newick_trees/science_hybrid_paper_species_phylogeny_TRIMMED_TIME_RENAMED_V2.tre")

#doesnt quite work as a fix as this doesnt have big birds
#auto <- auto1
#plot(auto)

auto<-read.tree("data/newick_trees/finches_324_shapeit4_phased_rh_maf0.15.recode_no_region_selection.for.poptree.txt.for.poptree.d3.b1.tre")

auto$tip.label <- gsub(".",". ",gsub("bigbirds","big bird lineage",auto$tip.label), fixed = T)

auto <- treeio::root(auto, c("L. noctis"), resolve.root = T, edgelabels = T)

auto$edge.length <- ifelse(auto$edge.length <0, 0, auto$edge.length)

mycalibration <- makeChronosCalib(auto, node="root", age.max=1)
mytimetree <- chronos(auto, lambda = 1, model = "correlated", calibration = mycalibration, control = chronos.control() )
write.tree(mytimetree, "output/newick_trees/auto_time_tree_V2.tre")
write.tree(as.newick(mytimetree), "output/newick_trees/auto_time_tree_V2_newick.tre")

finch.pruned.chrono <- read.tree(text=write.tree(mytimetree))
finch.pruned.chrono$tip.label <- gsub("_"," ", finch.pruned.chrono$tip.label)

finch.pruned.chrono3 <- phytools::reroot(finch.pruned.chrono, node = 24, position = .05)

library(phytools)
set.seed(998)

# output auto tree --------------------------------------------------------

ggtree(finch.pruned.chrono3) + geom_text(aes(label = node))

plot(rotateNodes(finch.pruned.chrono3,c(22,30,32,33,34,35,36,41,37,38,39,40)))

finch.pruned.chrono4 <- rotateNodes(finch.pruned.chrono3,c(22,30,32,33,34,35,36,41,37,38,39,40))

pdf("output/newick_trees/chrono_tree.pdf", width = 5.5, height = 8.5)
plot(finch.pruned.chrono4)
dev.off()
