library(tidyverse)
library(patchwork)
library(zoo)
library(ggrepel)
library(readxl)
library(GenomicRanges)
df.pixy <- read.csv("output/pixy/Four_geospiza_processed_Pixy.csv") %>% 
  select(-X)

#using filtered set of peak region
df.peaks.bed <- read.table("output/GEMMA/processed/autosomes_for_ful_mag_species_lmm_PEAKS.bed")


names(df.peaks.bed) <- c("chr","start","end","peak")
unique(df.peaks.bed$peak)

#add an extension around the peaks for plotting
df.peaks.bed <- df.peaks.bed %>% mutate(start.ext = start - 1000000, end.ext = end + 1000000)


df.pixy <- df.pixy %>% select(-contains("scandens"), -SNPs)

df.recomb <- read.csv("output/LDHelmet/windowed_rho_20kb.csv")

#with invariant sites
df.taj.d <- read.csv("data/tajimas_d/with_invariant/three_species_tajimasd_with_invariant.csv")

df.taj.d.long <- df.taj.d %>% select(-N_SNPS) %>% 
  pivot_longer(cols = -c("CHROM", "BIN_START", "BIN_END", "BIN_MID"), names_to = "stat", values_to = "value") %>% 
  mutate(Species = gsub("tajd_","", stat))



# format genes ------------------------------------------------------------
parv.df.an <- read.csv("data/annotation/Camarhynchus_parvulus_V1.1_genes.csv") 
df.genes.idx <- read.table("output/GEMMA/processed/autosomes_for_ful_mag_species_lmm_GENES.txt", header = T)
df.genes.idx <- df.genes.idx %>% select(gene_id, peak)
df.genes <- left_join(df.genes.idx, parv.df.an, by = "gene_id") %>% 
  dplyr::rename(CHROM = seqid) %>% mutate(mid = (start + end)/2) 

#convert gene symbols to uppercase (blast hits can be from species that are in lowercase on ensembl)
df.genes$gene_symbol <- toupper(df.genes$gene_symbol)

# load gemma data ---------------------------------------------------------

df.gemma.FOREGROUND <- read.table("output/GEMMA/processed/autosomes_for_ful_mag_species_lmm_FOR_PLOTTING_WITH_PIXY.txt", header = T)

# zfst --------------------------------------------------------------------

dfZ <- df.pixy %>%  filter(CHROM == "chrZ") %>% 
  mutate(zfst_fortis_fuliginosa = (fst_fortis_fuliginosa - mean(fst_fortis_fuliginosa, na.rm = T)) / sd(fst_fortis_fuliginosa, na.rm =T),
         zfst_fortis_magnirostris = (fst_fortis_magnirostris - mean(fst_fortis_magnirostris, na.rm = T)) / sd(fst_fortis_magnirostris, na.rm =T),
         zfst_fuliginosa_magnirostris = (fst_fuliginosa_magnirostris - mean(fst_fuliginosa_magnirostris, na.rm = T)) / sd(fst_fuliginosa_magnirostris, na.rm =T))

dfA <- df.pixy %>%  filter(CHROM != "chrZ") %>% 
  mutate(zfst_fortis_fuliginosa = (fst_fortis_fuliginosa - mean(fst_fortis_fuliginosa, na.rm = T)) / sd(fst_fortis_fuliginosa, na.rm =T),
         zfst_fortis_magnirostris = (fst_fortis_magnirostris - mean(fst_fortis_magnirostris, na.rm = T)) / sd(fst_fortis_magnirostris, na.rm =T),
         zfst_fuliginosa_magnirostris = (fst_fuliginosa_magnirostris - mean(fst_fuliginosa_magnirostris, na.rm = T)) / sd(fst_fuliginosa_magnirostris, na.rm =T))

df.overlayZFst <- rbind(dfA, dfZ)
df.pixy <- df.overlayZFst


# calc da -----------------------------------------------------------------


all_dxy <- grep("dxy_",names(df.pixy), value = T)

comp.dxy <- "dxy_fortis_fuliginosa"
#calculate da following Stankowski et al 2019 PLOS paper on monkeyflowers
for (comp.dxy in all_dxy){
  contr.spl <- str_split(comp.dxy, "_", simplify = T)
  pop1 <- paste(contr.spl[2], collapse="_") 
  pop2 <- paste(contr.spl[3], collapse="_") 
  df.pixy <- as.data.frame(df.pixy)
  df.pixy[,gsub("dxy","da",comp.dxy)] <- df.pixy[,comp.dxy] - rowMeans(df.pixy[,c(paste0("pi_",pop1),paste0("pi_",pop2))], na.rm = T)
  #df.pixy[,gsub("dxy","time",comp.dxy)] <- (df.pixy[,gsub("dxy","da",comp.dxy)] / (2*2.04e-9)) #calc with da
  df.pixy[,gsub("dxy","time",comp.dxy)] <- (df.pixy[,comp.dxy] / (2*2.04e-9)) #to calculate with dxy and not da
  #dont need to output mean pie directly
  #df.pixy[,gsub("dxy","meanpi",comp.dxy)] <- rowMeans(df.pixy[,c(paste0("pi_",pop1),paste0("pi_",pop2))], na.rm = T)
}


# plot without rollmean ---------------------------------------------------

df.pixy.long <- df.pixy %>% pivot_longer(cols = -c("CHROM","BIN_START","BIN_END","chr_ordered","row","chr_labels"), 
                                                             names_to = "stat", values_to = "value") %>% 
  mutate(contrast = gsub("da_","", gsub("time_","",gsub("fst_","", gsub("pi_","",gsub("zfst_","",gsub("dxy_","",stat))))))) %>% 
  mutate(BIN_MID = (BIN_START + BIN_END) / 2)

df4max <- df.pixy.long %>% filter(grepl("pi", stat))
df4max2 <- df.pixy.long %>% filter(grepl("zfst", stat))

#maxdxy <- max(df4max$value, na.rm = T)
maxdxy <- 0.006
maxzfst <- max(df4max2$value, na.rm =T)
maxzrecomb <- max(df.recomb$rho_per_kb, na.rm =T)
minzfst <- min(df4max2$value, na.rm =T)
maxlogp <- max(df.gemma.FOREGROUND$log_p, na.rm = T)

# get neutral background --------------------------------------------------
df.pixy.long <- df.pixy.long %>% mutate(index = paste(CHROM, BIN_START, BIN_END, sep="_"))

range.pixy <- GRanges(df.pixy.long$CHROM, IRanges(as.numeric(df.pixy.long$BIN_START), as.numeric(df.pixy.long$BIN_END)))
range.peaks_bed <- GRanges(df.peaks.bed$chr, IRanges(as.numeric(df.peaks.bed$start), as.numeric(df.peaks.bed$end)))
peaks.pixy.overlap <- as.data.frame(range.pixy[queryHits(findOverlaps(range.pixy, range.peaks_bed)), ])
peaks.pixy.overlap <- peaks.pixy.overlap %>% mutate(index = paste(seqnames, start, end, sep= "_"))
df.pixy.BACKGROUOND.long <- df.pixy.long %>% filter(!index %in% peaks.pixy.overlap$index) %>% 
  mutate(focus.pop = paste0(contrast, "_background"))

#for recombination 
df.recomb <- df.recomb %>% mutate(index = paste(CHROM, win_start, win_end, sep="_"))

range.recomb <- GRanges(df.recomb$CHROM, IRanges(as.numeric(df.recomb$win_start), as.numeric(df.recomb$win_end)))

peaks.recomb.overlap <- as.data.frame(range.recomb[queryHits(findOverlaps(range.recomb, range.peaks_bed)), ])
peaks.recomb.overlap <- peaks.recomb.overlap %>% mutate(index = paste(seqnames, start, end, sep= "_"))
df.recomb.BACKGROUOND.long <- df.recomb %>% filter(!index %in% peaks.recomb.overlap$index) %>% 
  mutate(focus.pop ="background")

#for tajimas d
df.taj.d.long <- df.taj.d.long %>% mutate(index = paste(CHROM, BIN_START, BIN_END, sep="_"))
range.taj.d <- GRanges(df.taj.d.long$CHROM, IRanges(as.numeric(df.taj.d.long$BIN_START), as.numeric(df.taj.d.long$BIN_END)))

peaks.taj.d.overlap <- as.data.frame(df.taj.d.long[queryHits(findOverlaps(range.taj.d, range.peaks_bed)), ])
peaks.taj.d.overlap <- peaks.taj.d.overlap %>%
  dplyr::mutate(index = paste(CHROM, BIN_START, BIN_END, sep= "_"))

df.taj.d.BACKGROUOND.long <- df.taj.d.long %>% filter(!index %in% peaks.recomb.overlap$index) %>% 
  mutate(focus.pop ="background")


# tables ------------------------------------------------------------------

#recombination

df.table.lists2 <- list()
for (peak.id in unique(df.peaks.bed$peak)){
  
  df.peak.focus <- df.peaks.bed %>% filter(peak == peak.id) %>% 
    mutate(mid = (start + end) /2)
  
  #create datafiles only for the region of interest for boxplots
  df.recomb.FOREGROUND <- df.recomb %>% 
    filter(CHROM == df.peak.focus$chr & win_start > df.peak.focus$start & win_end < df.peak.focus$end) %>% 
    mutate(focus.pop = "foreground")
  
  df.recomb.BACKGROUOND.sum <- df.recomb.BACKGROUOND.long %>%  
    summarise(median_value = median(rho_per_kb, na.rm = T))
  
  df4boxplot2 <- df.recomb.FOREGROUND %>% mutate(median_value = df.recomb.BACKGROUOND.sum$median_value)
  
  df.table.lists2[[peak.id]] <- df4boxplot2 %>% 
    mutate(value = rho_per_kb) %>% 
    mutate(delta_stat = value - median_value) %>% 
    summarise(median_stat = median(value, na.rm = T),
              mean_stat = mean(value, na.rm = T),
              sd_stat = sd(value, na.rm = T),
              n_stat = n(),
              median_delta_stat = median(delta_stat, na.rm = T),
              sd_delta_stat = sd(delta_stat, na.rm = T),
              .groups = "keep") %>%
    mutate(peak = peak.id, loc = paste0(df.peak.focus$chr, ":", round(df.peak.focus$start,0), "_", round(df.peak.focus$end,0)),
           se_stat = sd_stat / sqrt(n_stat),
           lower.ci_stat = mean_stat - qt(1 - (0.05 / 2), n_stat - 1) * se_stat,
           upper.ci_stat = mean_stat + qt(1 - (0.05 / 2), n_stat - 1) * se_stat,
           peak = peak.id) %>% 
    ungroup() 
}
all_peaks_stats2 = do.call(rbind, df.table.lists2)
mean(all_peaks_stats2$mean_stat)

df.table.lists <- list()
for (peak.id in unique(df.peaks.bed$peak)){
  
  df.peak.focus <- df.peaks.bed %>% filter(peak == peak.id) %>% 
    mutate(mid = (start + end) /2)
  
  df.rollmean.pixy2.FOREGROUND <- df.pixy.long %>% 
    filter(CHROM == df.peak.focus$chr & BIN_START > df.peak.focus$start & BIN_END < df.peak.focus$end) %>% 
    mutate(focus.pop = paste0(contrast, "_foreground"))
  
  ##The other approach is to calculate the change in stat vs the background
  df.pixy.BACKGROUOND.sum <- df.pixy.BACKGROUOND.long %>% group_by(stat) %>%  
    summarise(median_value = median(value, na.rm = T))
  
  df4boxplot <- df.rollmean.pixy2.FOREGROUND %>% left_join(df.pixy.BACKGROUOND.sum, by = "stat")
  
  #median value is the background
  df.table.lists[[peak.id]] <- df4boxplot %>% group_by(stat, contrast, median_value) %>% 
    mutate(delta_stat = value - median_value) %>% 
    summarise(median_stat = median(value, na.rm = T),
              mean_stat = mean(value, na.rm = T),
              sd_stat = sd(value, na.rm = T),
              n_stat = n(),
              median_delta_stat = median(delta_stat, na.rm = T),
              sd_delta_stat = sd(delta_stat, na.rm = T),
              .groups = "keep") %>%
    mutate(peak = peak.id, loc = paste0(df.peak.focus$chr, ":", round(df.peak.focus$start,0), "_", round(df.peak.focus$end,0)),
           se_stat = sd_stat / sqrt(n_stat),
           lower.ci_stat = mean_stat - qt(1 - (0.05 / 2), n_stat - 1) * se_stat,
           upper.ci_stat = mean_stat + qt(1 - (0.05 / 2), n_stat - 1) * se_stat,
           peak = peak.id,
           size = formatC((df.peak.focus$end - df.peak.focus$start), format = "d", big.mark = ",")) %>% 
    ungroup() %>% select(-median_value)
  
}

all_peaks_stats = do.call(rbind, df.table.lists)
#remove Z peaks, because foreground is determined by autosomal [i]
all_peaks_stats <- all_peaks_stats %>% filter(peak!=38 & peak!=39)


all_peaks_stats_wide <- all_peaks_stats %>% filter(!grepl("pi", stat) & stat != "no_sites") %>% 
  separate(stat, into = c("stat2",NA,NA), remove = T) %>%
  pivot_wider(names_from = stat2, values_from = c("median_stat", "sd_stat","median_delta_stat","sd_delta_stat",
                                                  "mean_stat", "se_stat", "lower.ci_stat", "upper.ci_stat")) 
  
all_peaks_stats_out <-  all_peaks_stats_wide %>% 
  select(peak, loc, size, contrast, median_stat_fst,median_stat_dxy,median_stat_da, median_stat_time, 
         fst = median_stat_fst, dxy = median_stat_dxy, da = median_stat_da, time = median_stat_time,
         region = loc) %>% 
  mutate(fst = round(fst, 2),
         dxy = round(dxy * 100, 2),
         da = round(da * 100, 2),
         time = round(time, 0))

name.updates <- read.csv("output/GEMMA/processed/old_new_peak_names.csv")
all_peaks_stats_out <- left_join(all_peaks_stats_out, name.updates, by = c("peak" = "peak.new")) 
#great.genes <- read.table("data/GREAT/Great_gene_by_region_simpl.bed") %>% select(V1, V4, peak = V1, gene = V4)
great.genes <- read.table("data/GREAT/GREAT_genes_by_regions_finches_20210702.txt") %>% select(V1, V2, peak = V2, gene = V1)

great.genes.w <- great.genes %>%
  group_by(peak) %>%
  summarise(Great.genes = toString(gene)) %>%
  ungroup() %>% mutate(contrast = "fortis_fuliginosa")
all_peaks_stats_out <- left_join(all_peaks_stats_out, great.genes.w, by = c("peak" = "peak", "contrast"))

df.all.genes <- df.genes %>% select(peak, gene_symbol) %>% 
  group_by(peak) %>% 
  summarise(all.gene_symbols = toString(gene_symbol)) %>% 
  mutate(contrast = "fortis_fuliginosa")

#how many (for text)
df.genes %>% 
  distinct(gene_symbol, .keep_all = T) %>% 
  group_by(peak) %>% summarise(n = n()) %>% print(n=28)
df.tmp <- df.genes %>% filter(peak == 26)
df.tmp %>% group_by(gene_symbol) %>% filter(n()>1)

df.all.genes.IDs <- df.genes %>% select(peak, gene_id) %>% 
  group_by(peak) %>% 
  summarise(all.gene_ensembl = toString(gene_id)) %>% 
  mutate(contrast = "fortis_fuliginosa")

all_peaks_stats_out <- left_join(all_peaks_stats_out, df.all.genes, by = c("peak" = "peak", "contrast"))
all_peaks_stats_out <- left_join(all_peaks_stats_out, df.all.genes.IDs, by = c("peak" = "peak", "contrast"))
all_peaks_stats_out <- all_peaks_stats_out %>% dplyr::rename(`Great genes` = Great.genes,
                                                      `Genes in region` = all.gene_symbols,
                                                      `Ensembl IDs` = all.gene_ensembl) %>% 
  select(-peak.original)

all_peaks_stats_out <- all_peaks_stats_out %>% select(-da) #dont output da
write_csv(all_peaks_stats_out, "assem_finch_code/supplementary_data_for_github/Supplementary_data_2_all_peaks_median_diversity_stat.csv", na = "")


#time
all_peaks_time <- all_peaks_stats %>% filter(grepl("time", stat)) %>% 
  mutate(loc = as.factor(loc), 
         peak_fact = as.factor(as.character(peak)))

all_peaks_time %>% filter(contrast == "fuliginosa_magnirostris") %>% 
  ggplot(aes(x = fct_reorder(peak_fact, -mean_stat), y = mean_stat/1000000)) +
  #geom_point(aes(x = loc, y = median_stat, fill = contrast, color = contrast), position = position_dodge(width=0.75)) +
  geom_point() +
  geom_linerange(aes(ymin=lower.ci_stat/1000000, ymax=upper.ci_stat/1000000), colour="black", width=.1, alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 10),
        legend.position = c(.8,.8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y = "Time (millions of years)", x = "Locus") +
  geom_hline(yintercept = .9, alpha = 0.5)
ggsave("output/pixy/all_peak_regions_time.pdf", width = 183, height = 95, useDingbats = F, units = "mm", dpi = "retina")
ggsave("output/pixy/all_peak_regions_time_sm.pdf", width = 90, height = 95, useDingbats = F, units = "mm", dpi = "retina")

#pi plots
all_peaks_pi <- all_peaks_stats %>% filter(grepl("pi", stat)) %>% mutate(species = contrast)
all_peaks_pi$species <- factor(all_peaks_pi$species, levels = c("fuliginosa", "fortis", "magnirostris"))
all_peaks_pi %>% 
  mutate(median_delta_stat = median_delta_stat * 1e3) %>% 
  ggplot() + 
  geom_boxplot(aes(x = species, y = median_delta_stat), outlier.size = -1) +
  geom_jitter(aes(x = species, y = median_delta_stat), width = 0.1, size = 0.9, alpha = 0.5, shape = 16) +
  theme_bw() + xlab(NULL) + 
  ylab(expression(paste(Delta, pi, " (x10"^3,")", sep = ""))) +
  theme(axis.text.x = element_text(face = "italic")) 
ggsave("output/pixy/all_peak_regions_delta_pi.pdf", width = 3, height = 3, useDingbats = F)
ggsave("output/pixy/all_peak_regions_delta_pi.png", width = 8, height = 6)

#how many peaks have the greatest reduction in diversity in this trio?
all_peaks_pi.wide <- all_peaks_pi %>% select(-species) %>% pivot_wider(id_cols = peak, names_from = "stat", values_from = "median_delta_stat")
all_peaks_pi.wide %>% filter(pi_magnirostris < pi_fortis & pi_magnirostris < pi_fuliginosa) %>% nrow() / nrow(all_peaks_pi.wide)

# per chr plots------------------------------------------------------------------------


chr.plots.pixy2 <- list()
df.lists <- list()
#V5 is after fixing log error
pdf(paste0("output/pixy/","three_species","_pixy_fst_dxy_pi_plot_no_rollmean_V5.pdf"), width = 8.5, height = 14, useDingbats = F)
for (peak.id in unique(df.peaks.bed$peak)){
  
  #create extended view for plotting
  
  df.peak.focus <- df.peaks.bed %>% filter(peak == peak.id)
  
  df.recomb.focus <- df.recomb %>% filter(CHROM == df.peak.focus$chr & win_start > df.peak.focus$start.ext & win_end < df.peak.focus$end.ext) %>% 
    mutate(focus.pop = "foreground")
  
  df.taj.d.focus <- df.taj.d.long %>% filter(CHROM == df.peak.focus$chr & BIN_START > df.peak.focus$start.ext & BIN_END < df.peak.focus$end.ext) 
  
  df.rollmean.pixy2.focus <- df.pixy.long %>% 
    filter(CHROM == df.peak.focus$chr & BIN_START > df.peak.focus$start.ext & BIN_END < df.peak.focus$end.ext) %>% 
    mutate(focus.pop = paste0(contrast, "_foreground"))
  
  df.gemma.FOREGROUND.focus <- df.gemma.FOREGROUND %>% 
    filter(chr == df.peak.focus$chr & ps > df.peak.focus$start.ext & ps < df.peak.focus$end.ext)
    
  
  zfst.rollmean.pixy2.focus <- df.rollmean.pixy2.focus %>% filter(grepl("zfst", stat))
  fst.rollmean.pixy2.focus <- df.rollmean.pixy2.focus %>% filter(grepl("fst", stat)) %>% filter(!grepl("zfst",stat))
  dxy.rollmean.pixy2.focus <- df.rollmean.pixy2.focus %>% filter(grepl("dxy", stat))
  pi.rollmean.pixy2.focus <- df.rollmean.pixy2.focus %>% filter(grepl("pi", stat))
  df.genes.focus <- df.genes %>% filter(peak == peak.id)
  
  general.theme <- theme(legend.position=c(0.9,0.8),
                         #panel.border=element_blank(),
                         panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line(),
                         #axis.title.x=element_blank(),
                         #axis.text.x = element_text(angle = 45, color = "black"),
                         #axis.text.x = element_blank(),
                         panel.grid = element_blank(),
                         panel.background = element_blank(),
                         panel.grid.major.y=element_blank(),
                         panel.grid.minor.y=element_blank(),
                         axis.title = element_text(size=14),
                         axis.text = element_text(size=14),
                         axis.ticks.x=element_blank(),
                         axis.ticks.y=element_line(size=0.2))
  
  #create datafiles only for the region of interest for boxplots
  df.recomb.FOREGROUND <- df.recomb %>% filter(CHROM == df.peak.focus$chr & win_start > df.peak.focus$start & win_end < df.peak.focus$end) %>% 
    mutate(focus.pop = "foreground")
  
  df.rollmean.pixy2.FOREGROUND <- df.pixy.long %>% 
    filter(CHROM == df.peak.focus$chr & BIN_START > df.peak.focus$start & BIN_END < df.peak.focus$end) %>% 
    mutate(focus.pop = paste0(contrast, "_foreground"))
  
  ###Two approaches to making these plots. One is to plot all foreground and background points
  #df4boxplot <- rbind(df.pixy.BACKGROUOND.long,df.rollmean.pixy2.FOREGROUND)
  #recomb4boxplot <- rbind(df.recomb.BACKGROUOND.long, df.recomb.FOREGROUND)
  
  ##The other approach is to calculate the change in stat vs the background
  df.pixy.BACKGROUOND.sum <- df.pixy.BACKGROUOND.long %>% group_by(stat) %>%  summarise(median_value = median(value, na.rm = T))
  df.recomb.BACKGROUOND.sum <- df.recomb.BACKGROUOND.long %>%  summarise(median_value = median(rho_per_kb, na.rm = T))
  df.taj.d.BACKGROUOND.long.sum <- df.taj.d.BACKGROUOND.long %>% group_by(stat) %>% summarise(median_value = median(value, na.rm = T))
  
  df4boxplot <- df.rollmean.pixy2.FOREGROUND %>% left_join(df.pixy.BACKGROUOND.sum, by = "stat")
  
  recomb4boxplot <- df.recomb.FOREGROUND %>% mutate(median_value = df.recomb.BACKGROUOND.sum$median_value)
  
  tajd4boxplot <- df.taj.d.focus %>% left_join(df.taj.d.BACKGROUOND.long.sum, by = "stat")
  
  df.pixy.long4zfstbox <- df4boxplot %>% filter(grepl("zfst", stat))
  df.pixy.long4dxybox <- df4boxplot %>% filter(grepl("dxy", stat))
  df.pixy.long4pibox <- df4boxplot %>% filter(grepl("pi", stat))
  
  #for the genes plot, if there are no genes it breaks the ggplot coommand
  #so if the genes df is empty, then just plot an empty rectangle
  if (dim(df.genes.focus)[1] == 0) {
    p.genes <- ggplot(data = zfst.rollmean.pixy2.focus, aes(x=BIN_MID/10000,y=value)) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y=element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank()) +
      ylim(1.9,2.6) +
      xlim(min(zfst.rollmean.pixy2.focus$BIN_MID, na.rm = T)/1000000, max(zfst.rollmean.pixy2.focus$BIN_MID,na.rm = T)/1000000)
  } else {
  
    #smart tip for aligning labels here: https://stackoverflow.com/questions/51024675/preserving-order-with-geom-text-repel
    #I try to add a bit of a buffer so that more can be plotted
    df.genes.focus$i <- seq(min(df.genes.focus$mid) - 500000, max(df.genes.focus$mid) + 500000, length.out = nrow(df.genes.focus))
    
    p.genes <- ggplot(data = zfst.rollmean.pixy2.focus, aes(x=BIN_MID/10000,y=value)) +
      geom_segment(data = df.genes.focus, aes(y = 2, yend = 2, x = (start - 50)/1000000, xend = (end + 50)/1000000), color = "black", size = 3) +
      geom_text(
        data = df.genes.focus,
        mapping = aes(y = 2.1, x = i/1000000, label = gene_symbol),
        parse = TRUE, hjust = 0, angle = 90, size = 3
      ) +
      geom_segment(
        data = df.genes.focus,
        mapping = aes(y = 2, yend = 2.1, x = mid/1000000, xend = i/1000000),
        size = 0.1
      ) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y=element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank()) +
      ylim(1.9,2.7) +
      xlim(min(zfst.rollmean.pixy2.focus$BIN_MID, na.rm = T)/1000000, max(zfst.rollmean.pixy2.focus$BIN_MID,na.rm = T)/1000000)
    
  #less effective with ggrepel
  #p.genes <- ggplot(data = zfst.rollmean.pixy2.focus, aes(x=BIN_MID/10000,y=value)) +
  #  geom_segment(data = df.genes.focus, aes(y = 2, yend = 2, x = (start - 50)/1000000, xend = (end + 50)/1000000), color = "black", size = 3) +
  #  geom_label_repel(data = df.genes.focus, 
  #                   aes(x = mid/1000000, y = 2, label = gene_symbol), max.overlaps = 20, nudge_y = .5, size = 2) +
  #  theme(panel.grid.major = element_blank(), 
  #        panel.grid.minor = element_blank(),
  #        panel.background = element_blank(), 
  #        axis.line = element_line(colour = "black"),
  #        axis.title.x=element_blank(),
  #        axis.ticks.x = element_blank(),
  #        axis.text.x = element_blank(),
  #        axis.title.y=element_blank(), 
  #        axis.text.y = element_blank(), 
  #        axis.ticks.y = element_blank()) +
  #  ylim(1.9,2.6) +
  #  xlim(min(zfst.rollmean.pixy2.focus$BIN_MID, na.rm = T)/1000000, max(zfst.rollmean.pixy2.focus$BIN_MID,na.rm = T)/1000000)
  #  #limits = c(min(df.subset$BIN_START, na.rm = T), max(df.subset$BIN_START, na.rm = T)
  }
  
  #points to highlight as outliers
  p.gemma.outlier <- df.gemma.FOREGROUND.focus %>%
    filter(chr == df.peak.focus$chr & ps > df.peak.focus$start & ps < df.peak.focus$end) %>% 
    filter(log_p > 7.3)
  ##First plot the gemma admixture result
  p.gemma <- ggplot() + 
    geom_point(data = df.gemma.FOREGROUND.focus, aes(x = ps/1000000, y = log_p), size = .5) +
    geom_point(data = p.gemma.outlier, aes(x = ps/1000000, y = log_p), size = .5, color = "red") +
    geom_vline(data = df.peak.focus, aes(xintercept = start/1000000), linetype = 2) +
    geom_vline(data = df.peak.focus, aes(xintercept = end/1000000), linetype = 2) +
    theme_bw() + 
    general.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank()) +
    labs(y=expression(paste("-log"[10], italic("P"),"-value")), x = unique(df.peak.focus$chr)) +
    #scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") + 
    ggtitle(paste0(unique(df.peak.focus$chr),", ", "peak = ",peak.id)) +
    ylim(0,maxlogp) 
  
  ##boxplot method 1 which plots both all background points and all foreground points
  #b.zfst <- ggplot() + 
  #  geom_boxplot(data = df.pixy.long4zfstbox, aes(x = focus.pop, y = value, fill = focus.pop), outlier.alpha = 0.1) +
  #  general.theme +
  #  theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank(),
  #        axis.title.y=element_blank(),
  #        legend.position = "none") +
  #  #labs(y=expression(ZF[ST])) +
  #  scale_fill_manual(values = c("grey80","#fc8d62","grey80","#66c2a5","grey80","#8da0cb"), name = "Contrast") 
  
  ##boxplot method2 which plots foreground points - median of all background points
  b.zfst <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_boxplot(data = df.pixy.long4zfstbox, aes(x = focus.pop, y = value - median_value, fill = focus.pop), outlier.alpha = 0.1) +
    general.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          #axis.title.y=element_blank(),
          legend.position = "none") +
    labs(y=expression(Delta)) +
    scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") 
  
  p.zfst <- ggplot() + 
    geom_line(data = zfst.rollmean.pixy2.focus, aes(x = BIN_MID/1000000, y = value, color = contrast), size = 1) +
    geom_vline(data = df.peak.focus, aes(xintercept = start/1000000), linetype = 2) +
    geom_vline(data = df.peak.focus, aes(xintercept = end/1000000), linetype = 2) +
    theme_bw() + 
    general.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank()) +
    labs(y=expression(italic(ZF[ST])), x = unique(df.peak.focus$chr)) +
    scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") + 
    ylim(minzfst,maxzfst) 
  
  p.zfst
  ##boxplot method 1 which plots both all background points and all foreground points
  #b.dxy <- ggplot() + 
  #  geom_boxplot(data = df.pixy.long4dxybox, aes(x = focus.pop, y = value, fill = focus.pop), outlier.alpha = 0.1) +
  #  general.theme +
  #  theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank(),
  #        axis.title.y=element_blank(),
  #        legend.position = "none") +
  #  #labs(y=expression(d[xy])) +
  #  scale_fill_manual(values = c("grey80","#fc8d62","grey80","#66c2a5","grey80","#8da0cb"), name = "Contrast")
  
  ##boxplot method2 which plots foreground points - median of all background points
  b.dxy <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_boxplot(data = df.pixy.long4dxybox, aes(x = focus.pop, y = value - median_value, fill = focus.pop), outlier.alpha = 0.1) +
    general.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          #axis.title.y=element_blank(),
          legend.position = "none") +
    labs(y=expression(Delta)) +
    scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") 
  
  #p.fst <- ggplot() + 
  #  geom_line(data = fst.rollmean.pixy2.focus, aes(x = BIN_MID/1000000, y = value, color = contrast), size = 1) +
  #  geom_vline(data = df.peak.focus, aes(xintercept = start/1000000), linetype = 2) +
  #  geom_vline(data = df.peak.focus, aes(xintercept = end/1000000), linetype = 2) +
  #  theme_bw() + 
  #  general.theme +
  #  theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank()) +
  #  labs(y=expression(F[ST]), x = unique(df.peak.focus$chr)) +
  #  scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") + 
  #  ggtitle(paste0("peak = ",peak.id)) +
  #  ylim(0,1) 
  
  p.dxy <- ggplot() + 
    geom_line(data = dxy.rollmean.pixy2.focus, aes(x = BIN_MID/1000000, y = value, color = contrast), size = 1) +
    geom_vline(data = df.peak.focus, aes(xintercept = start/1000000), linetype = 2) +
    geom_vline(data = df.peak.focus, aes(xintercept = end/1000000), linetype = 2) +
    theme_bw() + 
    general.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank()) +
    labs(y=expression(italic(d[XY])), x = unique(df.peak.focus$chr)) +
    scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") +
    ylim(0,maxdxy)
  
  #b.pi <- ggplot() + 
  #  geom_boxplot(data = df.pixy.long4pibox, aes(x = focus.pop, y = value, fill = focus.pop), outlier.alpha = 0.1) +
  #  general.theme +
  #  theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank(),
  #        axis.title.y=element_blank(),
  #        legend.position = "none") +
  #  #labs(y=expression(pi)) +
  #  scale_fill_manual(values = c("grey80","#fc8d62","grey80","#66c2a5","grey80","#8da0cb"), name = "Contrast")
  
  b.pi <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_boxplot(data = df.pixy.long4pibox, aes(x = focus.pop, y = value - median_value, fill = focus.pop), outlier.alpha = 0.1) +
    general.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          #axis.title.y=element_blank(),
          legend.position = "none") +
    labs(y=expression(Delta)) +
    scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") 
  
  
  p.pi <- ggplot() + 
    geom_line(data = pi.rollmean.pixy2.focus, aes(x = BIN_MID/1000000, y = value, color = contrast), size = 1) +
    geom_vline(data = df.peak.focus, aes(xintercept = start/1000000), linetype = 2) +
    geom_vline(data = df.peak.focus, aes(xintercept = end/1000000), linetype = 2) +
    theme_bw() + 
    general.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank()) +
    labs(y=expression(pi), x = paste0(unique(df.peak.focus$chr)," (Mb)")) +
    scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Population") +
    ylim(0,maxdxy)
  
  #b.rho <- ggplot() + 
  #  geom_boxplot(data = recomb4boxplot, aes(x = focus.pop, y = rho_per_kb, fill = focus.pop), outlier.alpha = 0.1) +
  #  general.theme +
  #  theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank(),
  #        axis.title.y=element_blank(),
  #        legend.position = "none") +
  #  labs(y="Rho/Kbp") +
  #  scale_fill_manual(values = c("grey80","orange"), name = "Contrast")
  
  b.rho <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_boxplot(data = recomb4boxplot, aes(x = focus.pop, y = rho_per_kb - median_value, fill = focus.pop), outlier.alpha = 0.1) +
    general.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          #axis.title.y=element_blank(),
          legend.position = "none") +
    labs(y=expression(Delta)) +
    scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") 
  
  p.rho <- ggplot() + 
    geom_line(data = df.recomb.focus, aes(x = win_mid/1000000, y = rho_per_kb, color = rho_per_kb), size = 1) +
    geom_vline(data = df.peak.focus, aes(xintercept = start/1000000), linetype = 2) +
    geom_vline(data = df.peak.focus, aes(xintercept = end/1000000), linetype = 2) +
    theme_bw() + 
    general.theme +
    #theme(legend.position=c(0.9,0.7)) +
    theme(legend.position="none") +
    labs(y=expression(paste(rho, "/Kb")), x = paste0(unique(df.peak.focus$chr)," (Mb)")) +
    scale_colour_gradient2(low="blue", high="red", mid = "orange", 
                           midpoint = (max(df.recomb$rho_per_kb, na.rm = TRUE) /2 ),
                           limits = c(0,maxzrecomb), name = "Rho/Kbp") +
    ylim(0,maxzrecomb)
  
  
  #### Tajimas
  b.tajd <- ggplot() + 
    geom_hline(yintercept = 0) +
    geom_boxplot(data = tajd4boxplot, aes(x = Species, y = value - median_value, fill = Species), outlier.alpha = 0.1) +
    general.theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          #axis.title.y=element_blank(),
          legend.position = "none") +
    labs(y=expression(Delta)) +
    scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") 
  
  p.tajd <- ggplot() + 
    geom_line(data = df.taj.d.focus, aes(x = BIN_MID/1000000, y = value, color = Species), size = 1) + 
    theme_bw() +
    general.theme +
    labs(y="Tajima's D", x = paste0(unique(df.peak.focus$chr)," (Mb)")) +
    geom_vline(data = df.peak.focus, aes(xintercept = start/1000000), linetype = 2) +
    geom_vline(data = df.peak.focus, aes(xintercept = end/1000000), linetype = 2) +
    scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Population")
    
  
  layout <- "
  IIIIIII##
  JJJJJJJ##
  JJJJJJJ##
  AAAAAAABB
  AAAAAAABB
  AAAAAAABB
  CCCCCCCDD
  CCCCCCCDD
  CCCCCCCDD
  EEEEEEEFF
  EEEEEEEFF
  EEEEEEEFF
  GGGGGGGHH
  GGGGGGGHH
  "
  
  #plot with recombination. Move GGGGGGGHH to one line
  #  chr.plots.pixy2[[peak.id]] <- p.zfst + b.zfst + p.dxy +b.dxy + p.pi + b.pi +  p.rho + b.rho +p.genes + p.gemma +
  #    plot_layout(design = layout)
  
  #plot with tajimas D instead
  chr.plots.pixy2[[peak.id]] <- p.zfst + b.zfst + p.dxy +b.dxy + p.pi + b.pi +  p.tajd + b.tajd +p.genes + p.gemma +
    plot_layout(design = layout)

  print(chr.plots.pixy2[[peak.id]])
  
  #output a dataframe of pi 
  df.lists[[peak.id]] <- df.pixy.long4pibox %>% group_by(stat, median_value) %>% 
    summarise(median_pi = median(value, na.rm = T),
              .groups = "keep") %>%
    mutate(delta_pi = median_pi - median_value) %>% 
    mutate(peak = peak.id)
}
dev.off()

pdf("output/pixy/peaks_3_for_supplement.pdf", width = 8.5, height = 14, useDingbats = F)
(chr.plots.pixy2[[3]])
dev.off()


pdf("output/pixy/peaks_26_for_supplement.pdf", width = 8.5, height = 14, useDingbats = F)
(chr.plots.pixy2[[26]])
dev.off()

pdf("output/pixy/peaks_24_for_supplement.pdf", width = 8.5, height = 14, useDingbats = F)
(chr.plots.pixy2[[24]])
dev.off()

# investigate chr1A -------------------------------------------------------
df.gemma.chr1A <- read.table("output/GEMMA/processed/chr1A_for_ful_mag_species_lmm_FOR_PLOTTING_WITH_PIXY.txt", header = T)

p1 <- ggplot(df.gemma.chr1A) + geom_point(aes(x = ps, y = log_p),size = 1) + theme_bw()
p2 <- df.recomb %>% filter(CHROM == "chr1A") %>% 
  ggplot() + geom_line(aes(x = win_mid/1000000, y = rho_per_kb)) + 
  theme_bw() +
  labs(y=expression(paste(rho, "/Kb"), sep = ""), x = NULL) +
  theme(text = element_text(size = 18))
  #xlim(30, 45)

#p2
p3 <- df.pixy.long %>% 
  filter(CHROM == "chr1A" & grepl("zfst",stat)) %>% 
  ggplot() + 
  geom_line(aes(x = BIN_MID/1000000, y = value, color = contrast)) + 
  theme_bw() + 
  scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb")) +
  theme(legend.position = "none") +
  labs(y=expression(italic(ZF[ST])), x = "chr1A (Mb)") +
  theme(text = element_text(size = 18))
  #xlim(30, 45)


p2/p3

#p1/p2/p3
#ggsave("output/pixy/chr1A_compare_recombination_rate_fst_pvalue.png", width = 12, height = 8)
ggsave("output/pixy/chr1A_compare_recombination_rate_fst_pvalue_V2.png", width = 12, height = 4)
