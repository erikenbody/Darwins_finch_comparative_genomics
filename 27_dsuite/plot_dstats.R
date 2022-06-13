library(tidyverse)


df.dstats <- read.table("data/dsuite/all_finch_trios_BBAA.txt", header = T)
df.prop.mag <- df.dstats %>% filter(P3 == "Genovesa_magnirostris" & P2 == "Genovesa_propinqua")
df.mag.prop <- df.dstats %>% filter(P2 == "Genovesa_magnirostris" & P3 == "Genovesa_propinqua")


# fd ----------------------------------------------------------------------


df.fd.prop.mag.pinta <- read.table("data/dsuite/Pinta_scandens_Genovesa_propinqua_Genovesa_magnirostris_localFstats_prop_and_mag2_5025_50_25.txt", header = T)
#df.fd.prop.mag.pinta %>% 
#  filter(chr == "chr9") %>% 
#  mutate(f_d = ifelse(f_d < 0, 0, f_d)) %>% 
#  ggplot() + geom_line(aes(x = windowStart/1000000, y = f_d)) +
#  xlim(9,12) +
#  geom_vline(xintercept = 10002322/1000000)
#
#df.fd.prop.mag.pinta %>% 
#  filter(chr == "chr1") %>% 
#  mutate(f_d = ifelse(f_d < 0, 0, f_d)) %>% 
#  ggplot() + geom_point(aes(x = windowStart/1000000, y = d_f)) +  
#  geom_vline(xintercept = 62446429/1000000, color = "red") +xlim(60,65)
#
#df.fd.prop.mag.daphne <- read.table("data/dsuite/Daphne_scandens_Genovesa_propinqua_Genovesa_magnirostris_localFstats_prop_and_mag2_100_50.txt", header = T)
#df.fd.prop.mag.daphne %>% filter(chr == "chr9") %>% 
#  #mutate(f_d = ifelse(f_d < 0, 0, f_d)) %>% 
#  ggplot() + geom_line(aes(x = windowStart/1000000, y = f_d)) +
#  xlim(10,20)
#
#
#df.fd.prop.mag.daphne %>% 
#  filter(chr == "chr2") %>% 
#  #mutate(f_d = ifelse(f_d < 0, 0, f_d)) %>% 
#  ggplot() + geom_line(aes(x = windowStart/1000000, y = f_d)) +
#  geom_vline(xintercept = 62446429/1000000, color = "red") + ylim(0,3)
#

# summarize ---------------------------------------------------------------

library(tidyverse)
library(patchwork)


#using filtered set of peak region
df.peaks.bed <- read.table("output/GEMMA/processed/autosomes_for_ful_mag_species_lmm_PEAKS.bed")


names(df.peaks.bed) <- c("chr","start","end","peak")
unique(df.peaks.bed$peak)

#add an extension around the peaks for plotting
df.peaks.bed <- df.peaks.bed %>% mutate(start.ext = start - 1000000, end.ext = end + 1000000)

# get neutral background --------------------------------------------------
df.fd.prop.mag.pinta.long <- df.fd.prop.mag.pinta %>% pivot_longer(cols = -c("chr","windowStart", "windowEnd"), names_to = "stat")


df.fd.prop.mag.pinta.long <- df.fd.prop.mag.pinta.long %>% 
  mutate(index = paste(chr, windowStart, windowEnd, sep="_")) %>% 
  mutate(contrast = "prop.mag.pinta")


range.pixy <- GRanges(df.fd.prop.mag.pinta.long$chr, IRanges(as.numeric(df.fd.prop.mag.pinta.long$windowStart), as.numeric(df.fd.prop.mag.pinta.long$windowEnd)))
range.peaks_bed <- GRanges(df.peaks.bed$chr, IRanges(as.numeric(df.peaks.bed$start), as.numeric(df.peaks.bed$end)))
peaks.pixy.overlap <- as.data.frame(range.pixy[queryHits(findOverlaps(range.pixy, range.peaks_bed)), ])
peaks.pixy.overlap <- peaks.pixy.overlap %>% mutate(index = paste(seqnames, start, end, sep= "_"))
df.pixy.BACKGROUOND.long <- df.fd.prop.mag.pinta.long %>% filter(!index %in% peaks.pixy.overlap$index) %>%
  mutate(contrast = "prop.mag.pinta") %>% 
  mutate(focus.pop = paste0(contrast, "_background"))

##The other approach is to calculate the change in stat vs the background
df.pixy.BACKGROUOND.sum <- df.pixy.BACKGROUOND.long %>% group_by(stat) %>%  
  summarise(median_value = median(value, na.rm = T),
            mean_stat = mean(value, na.rm = T),
            n_stat = n(),
            sd_stat = sd(value, na.rm = T)) %>% 
  mutate(se_stat = sd_stat / sqrt(n_stat),
         lower.ci_stat = mean_stat - qt(1 - (0.05 / 2), n_stat - 1) * se_stat,
         upper.ci_stat = mean_stat + qt(1 - (0.05 / 2), n_stat - 1) * se_stat)


# tables ------------------------------------------------------------------



df.table.lists <- list()
for (peak.id in unique(df.peaks.bed$peak)){
  
  df.peak.focus <- df.peaks.bed %>% filter(peak == peak.id) %>% 
    mutate(mid = (start + end) /2)
  
  df.rollmean.pixy2.FOREGROUND <- df.fd.prop.mag.pinta.long %>% 
    filter(chr == df.peak.focus$chr & windowStart > df.peak.focus$start & windowEnd < df.peak.focus$end) %>% 
    mutate(focus.pop = paste0(contrast, "_foreground"))
  
  
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

all_peaks_stats_wide <- all_peaks_stats %>% filter(!grepl("pi", stat) & stat != "no_sites") %>% 
  pivot_wider(names_from = stat, values_from = c("median_stat", "sd_stat","median_delta_stat","sd_delta_stat",
                                                  "mean_stat", "se_stat", "lower.ci_stat", "upper.ci_stat")) 

all_peaks_stats_wide %>% ggplot() + geom_point(aes(x = as.factor(peak), y = median_delta_stat_d_f))

#all_peaks_stats_out <-  all_peaks_stats_wide %>% 
#  select(peak, loc, size, contrast, median_stat_fst,median_stat_dxy,median_stat_da, median_stat_time, 
#         fst = median_stat_fst, dxy = median_stat_dxy, da = median_stat_da, time = median_stat_time,
#         region = loc) %>% 
#  mutate(fst = round(fst, 2),
#         dxy = round(dxy * 100, 2),
#         da = round(da * 100, 2),
#         time = round(time, 0))
#

p.intro <- all_peaks_stats_wide %>% 
  mutate(peak_fact = as.factor(peak)) %>% 
  #ggplot(aes(x = fct_reorder(peak_fact, -mean_stat_d_f), y = mean_stat_d_f)) +#plot in order of highest value
  ggplot(aes(x = peak_fact, y = mean_stat_d_f)) +#plot in order of highest value
  #geom_point(aes(x = loc, y = median_stat, fill = contrast, color = contrast), position = position_dodge(width=0.75)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower.ci_stat_d_f, ymax=upper.ci_stat_d_f), colour="black", width=.1, alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 10),
        legend.position = c(.8,.8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y = expression(paste('Fraction of introgression ', d[f])), x = "Locus") +
  geom_hline(yintercept = 0.0435, alpha = 0.5) #+ #median d_f across background regions
  #geom_rect(aes(xmin = 0, xmax = 28, ymin = 0.0427, ymax = 0.0442), color = "red")

#ggsave("output/dstats/all_peak_regions_propinqua_magnirostris_introgression.pdf", width = 90, height = 95, useDingbats = F, units = "mm", dpi = "retina")
#p1 is from allele frequency script
p1 / p.intro
ggsave("output/dstats/all_peak_regions_propinqua_magnirostris_introgression_with_af.pdf", width = 12, height = 8)

#small for main plot
all_peaks_stats_wide %>% 
  mutate(peak_fact = as.factor(peak)) %>% 
  #ggplot(aes(x = fct_reorder(peak_fact, -mean_stat_d_f), y = mean_stat_d_f)) +#plot in order of highest value
  ggplot(aes(x = fct_reorder(peak_fact, -mean_stat_d_f), y = mean_stat_d_f)) +#plot in order of highest value
  #geom_point(aes(x = loc, y = median_stat, fill = contrast, color = contrast), position = position_dodge(width=0.75)) +
  geom_point() +
  geom_linerange(aes(ymin=lower.ci_stat_d_f, ymax=upper.ci_stat_d_f), colour="black", width=.1, alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 10),
        legend.position = c(.8,.8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y = expression(paste('Fraction of introgression ', d[f])), x = "Locus") +
  geom_hline(yintercept = 0.0435, alpha = 0.5) #+ #median d_f across background regions
ggsave("output/dstats/all_peak_regions_df_sm.pdf", width = 90, height = 95, useDingbats = F, units = "mm", dpi = "retina")



# man ---------------------------------------------------------------------


source("assem_finch_code/07_popgen/generic_manhattan_fnc.R")
#changed to log10 june 6
fd.order <- function(gemma.path, input.var){
  df.gemma <- gemma.path
  df.gemma$chr_ordered <- factor(df.gemma$chr, levels = chr_order)
  df.gemma <- df.gemma %>% dplyr::arrange(chr_ordered, windowStart)
  df.gemma$row<-1:nrow(df.gemma)
  df.gemma <- as.data.frame(df.gemma)
  #df.gemma$log_p <- -log10(df.gemma[,input.var])
  df.gemma$introgression_rollmean <- zoo::rollmean(df.gemma[,input.var],50,fill=NA)
  df.gemma$chr_labels <- gsub("chr", "", df.gemma$chr_ordered)
  chr_breaks <- df.gemma %>% group_by(chr_labels) %>% dplyr::summarise(chr_breaks = mean(row))
  df.gemma
}

manc3 <- function(df.in, input.var){
  chr_breaks <- df.in %>% filter(chr_ordered != "chrunknown" & !is.na(row)) %>% 
    mutate(chr_ordered = factor(chr_ordered, levels = chr_order)) %>%
    group_by(chr_ordered, chr_labels) %>% 
    dplyr::summarise(chr_breaks = mean(row))
  
  chrom.colors <- data.frame(chr_ordered = grep("chr", unique(df.in$chr_ordered), value = T),
                             color.num = rep(1:2,length(grep("chr", unique(df.in$chr_ordered))))) %>% 
    distinct(chr_ordered, .keep_all = T)
  
  df.in2 <- df.in %>% #mutate(row = 1:n()) %>% 
    left_join(chrom.colors, by = "chr_ordered") %>% 
    mutate(color.num = as.factor(color.num))
  
  df.in2 %>% filter(chr_labels != "unknown" & !is.na(row)) %>%
    ggplot(aes_string(x = "row", y = input.var, col = "color.num")) + theme_bw() +
    theme(legend.position="none",
          #panel.border=element_blank(),
          panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line(),
          axis.title.x=element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(size=10),
          axis.text.x = element_text(size = 6.5),
          axis.text = element_text(size=10),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    geom_point(size=0.9,shape=20,stroke=0.2) +
    scale_color_manual(values=rep(c("grey30","grey70"))) +
    #scale_color_manual(values = c(rep_len(c("grey30", "red"), length(unique(chr_breaks$chr_ordered))+1))) #+
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                       labels = function(labels) {
                         sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                       }) +
    labs(y = expression(paste('Fraction of introgression ', d[f])), x = NULL) 
  #scale_x_continuous(breaks=chr_breaks$chr_breaks, 
  #                   labels = chr_breaks$chr_labels)
}


df.fd.prop.mag.pinta.order <- fd.order(df.fd.prop.mag.pinta, "d_f")
df.fd.prop.mag.pinta.order <- df.fd.prop.mag.pinta.order %>%  mutate(index = paste(chr, windowStart, windowEnd, sep="_")) %>% 
  mutate(peak = ifelse(!index %in% df.pixy.BACKGROUOND.long$index, "foreground", "background"))

df.fd.prop.mag.pinta.outliers <- df.fd.prop.mag.pinta.order %>% filter(peak == "foreground")

png("output/dstats/whole_genome_prop.mag.pinta.png", width = 2400, height = 400, res = 200) 
manc3(df.fd.prop.mag.pinta.order, "d_f") +
  geom_point(data = df.fd.prop.mag.pinta.outliers, aes(x = row, y = d_f), color = "red", size=0.9,shape=20,stroke=0.2)
dev.off()

p.man.fd <- manc3(df.fd.prop.mag.pinta.order, "d_f") +
  geom_point(data = df.fd.prop.mag.pinta.outliers, aes(x = row, y = d_f), color = "red", size=0.9,shape=20,stroke=0.2)

#have to run af script and the propinqua dxy script
af.plot <- p1 +  theme(legend.position = "bottom")
pdf("output/dstats/all_peak_regions_propinqua_magnirostris_introgression_with_af_man.pdf",width = 8.5, height = 11)
p.man.fd /  p.intro / p.prop.dxy / af.plot
dev.off()


