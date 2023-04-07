library(GenomicRanges)
library(patchwork)
library(gghighlight)
library(tidyverse)

source("assem_finch_code/17_pixy/pixy_manhattan_fnc.R")

#pixy_contrast <- "data/pixy/fuliginosa_fortis/output/"

#df.tmp <- read.table("data/pixy/four_species/pixy_four_geospiza_6_fst.txt", header = T)
#df.tmp %>% filter(window_pos_1 == 24780001)

df.pixy_out <- input_pixy_V2("data/pixy/four_species/")

#70 windows with < 100 SNPs
df.pixy_out <- df.pixy_out %>% filter(SNPs > 100) 

write.csv(df.pixy_out,"output/pixy/Four_geospiza_processed_Pixy.csv")

all_contrasts <- gsub("fst_","",grep("fst_",names(df.pixy_out), value = T))

#make a pdf for each contrast per chromosome
for (contrast_name in all_contrasts){
  print(contrast_name)
  plot_pixy_V2(df.pixy_out, contrast_name)
}

# load and process 1kbp windows -------------------------------------------


#df.pixy_out_1kb <- input_pixy_V2("data/pixy/four_species_1kb/")
#write.csv(df.pixy_out_1kb,"output/pixy/Four_geospiza_processed_Pixy_1kb.csv")
#
#df.pixy_out_1kb.pi <- df.pixy_out_1kb %>% select(CHROM, BIN_START, BIN_END, contains("pi"))
#write.csv(df.pixy_out_1kb.pi,"output/pixy/Four_geospiza_processed_Pixy_1kb_pi.csv")

# manhattan plot ----------------------------------------------------------

p1 <- manc(df.pixy_out, "fst_fortis_fuliginosa") + ggtitle("fst_fortis_fuliginosa")
p2 <- manc(df.pixy_out, "fst_fortis_magnirostris") + ggtitle("fst_fortis_magnirostris")
p3 <- manc(df.pixy_out, "fst_fortis_scandens") + ggtitle("fst_fortis_scandens")
p4 <- manc(df.pixy_out, "fst_fuliginosa_magnirostris") + ggtitle("fst_fuliginosa_magnirostris")
p5 <- manc(df.pixy_out, "fst_fuliginosa_scandens") + ggtitle("fst_fuliginosa_scandens")
p6 <- manc(df.pixy_out, "fst_magnirostris_scandens") + ggtitle("fst_magnirostris_scandens")

p1 / p2 / p3 / p4 / p5 / p6
ggsave("output/pixy/six_fst_plots.png", height = 11, width = 8.5)


# overlay -----------------------------------------------------------------

df.overlay <- df.pixy_out %>% select(CHROM, BIN_START, BIN_END, chr_ordered, row, chr_labels, 
                       fst_fortis_fuliginosa, fst_fortis_magnirostris, fst_fuliginosa_magnirostris)

dfZ <- df.overlay %>%  filter(CHROM == "chrZ") %>% 
  mutate(zfst_fortis_fuliginosa = (fst_fortis_fuliginosa - mean(fst_fortis_fuliginosa, na.rm = T)) / sd(fst_fortis_fuliginosa, na.rm =T),
         zfst_fortis_magnirostris = (fst_fortis_magnirostris - mean(fst_fortis_magnirostris, na.rm = T)) / sd(fst_fortis_magnirostris, na.rm =T),
         zfst_fuliginosa_magnirostris = (fst_fuliginosa_magnirostris - mean(fst_fuliginosa_magnirostris, na.rm = T)) / sd(fst_fuliginosa_magnirostris, na.rm =T))

dfA <- df.overlay %>%  filter(CHROM != "chrZ") %>% 
  mutate(zfst_fortis_fuliginosa = (fst_fortis_fuliginosa - mean(fst_fortis_fuliginosa, na.rm = T)) / sd(fst_fortis_fuliginosa, na.rm =T),
         zfst_fortis_magnirostris = (fst_fortis_magnirostris - mean(fst_fortis_magnirostris, na.rm = T)) / sd(fst_fortis_magnirostris, na.rm =T),
         zfst_fuliginosa_magnirostris = (fst_fuliginosa_magnirostris - mean(fst_fuliginosa_magnirostris, na.rm = T)) / sd(fst_fuliginosa_magnirostris, na.rm =T))

df.overlayZFst <- rbind(dfA, dfZ)
df.overlayZFst <- df.overlayZFst %>% select(CHROM, BIN_START, BIN_END, chr_ordered, row, chr_labels, 
                          zfst_fortis_fuliginosa, zfst_fortis_magnirostris, zfst_fuliginosa_magnirostris)

df.overlay.long <- df.overlayZFst %>% pivot_longer(cols = -c("CHROM","BIN_START","BIN_END","chr_ordered","row","chr_labels"), 
                            names_to = "contrast", values_to = "Zfst")

df.overlay.long <- df.overlay.long %>% 
  mutate(contrast_label = ifelse(contrast == "zfst_fortis_fuliginosa", "fortis vs. fuliginosoa",
                                 ifelse(contrast == "zfst_fortis_magnirostris", "fortis vs. magnirostris",
                                        ifelse(contrast == "zfst_fuliginosa_magnirostris", "fuliginosa vs. magnirostris", NA))))
                           
chr_breaks <- df.overlay.long %>% filter(chr_ordered != "chrunknown" & !is.na(row)) %>% 
  mutate(chr_ordered = factor(chr_ordered, levels = chr_order)) %>%
  group_by(chr_ordered, chr_labels) %>% 
  dplyr::summarise(chr_breaks = mean(row), 
                   chr_ends = max(row),.groups = "keep")

p.overlay <- df.overlay.long %>% 
  #filter(CHROM!="chrZ") %>% #for now remove Z because not plotted with GEMMA
  ggplot() + 
  geom_point(aes(x = row, y = Zfst, color = contrast_label), alpha = 0.5, size=1.5,shape=20,stroke=0.2) + 
  theme_bw() + 
  theme(legend.position=c(0.9,0.8),
        #panel.border=element_blank(),
        panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line(),
        axis.title.x=element_blank(),
        #axis.text.x = element_text(angle = 45, color = "black"),
        #axis.text.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        axis.title.y = element_text(size=10),
        axis.text = element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2))+
        #legend.background = element_rect(fill=NA)) +
  labs(y=expression(ZF[ST])) +
  scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") +
  scale_alpha(guide = 'none') +
  scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                     labels = function(labels) {
                       sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                     }) +
  labs(y=expression(ZF[ST])) +
  geom_vline(data = chr_breaks, aes(xintercept = chr_ends), alpha = 0.1) +
  guides(colour = guide_legend(override.aes = list(size=10)),
         shape = guide_legend(override.aes = list(size = 10)))

p.overlay
ggsave("output/pixy/three_species_overlay_Zfst.png", width = 12, height = 4, dpi = 300)
p.overlay + scale_y_reverse() + theme(legend.position=c(0.9,0.3))
ggsave("output/pixy/three_species_overlay_Zfst_reverse.png", width = 12, height = 4, dpi = 300)


# lines? ------------------------------------------------------------------

drop.these.chr.labs <- c("chr21","chr22","chr23","chr24","chr25","chr26", "chr27","chr28", "chr29")
chr_breaks <- chr_breaks %>% mutate(chr_labels = ifelse(chr_ordered %in% drop.these.chr.labs, "", chr_labels))
p.overlay2 <- df.overlay.long %>% 
  mutate(contrast_label = gsub("vs.","-", contrast_label)) %>%
  #filter(CHROM!="chrZ") %>% #for now remove Z because not plotted with GEMMA
  ggplot() + 
  geom_line(aes(x = row, y = Zfst, color = contrast_label), alpha = 0.5, size=.5) + 
  theme_bw() + 
  theme(legend.position=c(0.9,0.8),
        #panel.border=element_blank(),
        panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line(),
        axis.title.x=element_blank(),
        #axis.text.x = element_text(angle = 45, color = "black"),
        #axis.text.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        axis.title.y = element_text(size=6),
        axis.text = element_text(size=6),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2),
        legend.text = element_text(face = "italic", size = 6),
        legend.title = element_text(size = 6),
        legend.margin=margin(t=-10))+
  #legend.background = element_rect(fill=NA)) +
  labs(y=expression(ZF[ST])) +
  scale_color_manual(values = c("#fc8d62","#66c2a5","#8da0cb"), name = "Contrast") +
  scale_alpha(guide = 'none') +
  scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                     labels = function(labels) {
                       sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                     }) +
  labs(y=expression(italic(ZF[ST]))) +
  geom_vline(data = chr_breaks, aes(xintercept = chr_ends), alpha = 0.1) +
  guides(colour = guide_legend(override.aes = list(size=5)),
         shape = guide_legend(override.aes = list(5)))


p.overlay2
ggsave("output/pixy/three_species_overlay_Zfst_lines.png", width = 12, height = 4, dpi = 300)
p.overlay2 + scale_y_reverse() + theme(legend.position="bottom")
ggsave("output/pixy/three_species_overlay_Zfst_reverse_lines.pdf", width = 183, height = 61, dpi = 300, units = "mm")


