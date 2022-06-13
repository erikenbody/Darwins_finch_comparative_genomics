library(tidyverse)
library(patchwork)
library(readxl)
df.af <- read.csv("data/allele_frequencies/GenotypeDiffs_28_loci.csv")
df.islm <- read_excel("data/morphology/Geospiza + tree finches.xlsx")
df.peaks.bed <- read.table("output/GEMMA/processed/autosomes_for_ful_mag_species_lmm_PEAKS.bed")

df.af.long <- df.af %>% 
  pivot_longer(-c("Species", "Order"), names_to = "Locus", values_to = "AF") %>% 
  mutate(Locus = as.factor(as.numeric(gsub("Locus","",Locus))))

p1 <- ggplot() + 
  geom_tile(data = df.af.long, aes(x = Locus, y = fct_reorder(Species, -Order), fill = AF)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(face = 'italic'),
        legend.position = "top") +
  labs(y  = NULL) +
  scale_fill_gradient2(low="steelblue2", high="tomato3", mid = "white", midpoint = (max(df.af.long$AF, na.rm = TRUE) /2 ), name = "DAF") 
p1
ggsave("output/allele_frequency/heatmap_daf.pdf", width = 8, height = 6)  


ggplot(data = df.af.long, aes(x = Locus, y = fct_reorder(Species, -Order), fill = AF)) + 
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(face = 'italic'),
        legend.position = "top") +
  labs(y  = "Species") +
  geom_text(aes(label = round(AF, 1)), size = 2) +
  scale_fill_gradient2(low="steelblue2", high="tomato3", mid = "white", midpoint = (max(df.af.long$AF, na.rm = TRUE) /2 ), name = "DAF") 
ggsave("output/allele_frequency/heatmap_daf_text.pdf", width = 8, height = 6)  


# morphology --------------------------------------------------------------
df.order <- df.af %>% select(Species, Order)
df.islm <- left_join(df.islm, df.order, by = c("species" = "Species")) %>% 
  filter(!is.na(Order))

df.islm.sum <- df.islm %>% group_by(species, Order) %>% 
  summarise(mean_pca = mean(`PC1 beak`, na.rm = T)) %>% 
  mutate(morphology = "Bill size")


p2 <- ggplot(data = df.islm.sum, aes(x = morphology, y = fct_reorder(species, -Order), fill = mean_pca)) + 
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0),
        axis.text.y = element_blank(),
        legend.position = "top") +
  scale_fill_gradient2(high="steelblue2", low="tomato3", mid = "white", na.value = "grey80", midpoint = (max(df.af.long$AF, na.rm = TRUE) /2 ), name = "PC1")  +
  ylab(NULL) + xlab(NULL)

p1 + p2 + plot_layout(widths = c(10,1)) +
  plot_layout(guides = "collect") & theme(legend.position = 'top')
ggsave("output/allele_frequency/heatmap_daf_text_morphology.pdf", width = 8, height = 6)  

##length
df.peaks.bed <- df.peaks.bed %>% mutate(length = V3 - V2) %>% 
  dplyr::rename(peak = V4) %>% 
  mutate(chr.id = as.factor(gsub("chr","",V1)))

p3 <- ggplot(df.peaks.bed) + 
  geom_point(aes(x = as.factor(peak), y = length/1000000)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        text = element_text(size = 10),
        legend.position = c(.8,.8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y = NULL, x = NULL)

layout <- "
AAAAAAAAAA#
BBBBBBBBBBC
BBBBBBBBBBC
BBBBBBBBBBC
BBBBBBBBBBC
BBBBBBBBBBC
"


p3 + p1 + p2 + 
  plot_layout(design = layout, guides = "collect") & theme(legend.position = "bottom")
ggsave("output/allele_frequency/heatmap_daf_text_morphology_V2.pdf", width = 8, height = 6)  


p3 + p1 + p2 + 
  plot_layout(design = layout, guides = "collect") & theme(legend.position = "none")
ggsave("output/allele_frequency/heatmap_daf_text_morphology_no_legend.pdf", width = 8, height = 6)  


