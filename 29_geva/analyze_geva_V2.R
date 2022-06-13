library(tidyverse)
library(data.table)
library(patchwork)
library(ggridges)


df.top100<-read_tsv("output/GEMMA/processed/autosomes_for_ful_mag_species_lmm_OUTLIER_SNPs_TOP100_truly100.txt",
                    col_names = F)

for (peak in unique(df.top100$X4)){
  df.peak <- df.top100 %>% 
    filter(X4 == peak) %>% 
    select(X2) %>% 
    write_tsv(paste0("data/geva/",peak, "_positions.txt"), col_names = F)
}

##make bed file with extensions
df.bed <- read.table("output/GEMMA/processed/autosomes_for_ful_mag_species_lmm_PEAKS.bed")
df.bed$start <- df.bed$V2 - 50000
df.bed$end <- df.bed$V3 + 50000

df.bed %>% 
  select(V1, start, end, V4) %>% 
  write_tsv("output/GEMMA/processed/autosomes_for_ful_mag_species_lmm_PEAKS_extension.bed", col_names = F)

# load into one datafile --------------------------------------------------

file.path.geva <- grep("*sites.txt",list.files("data/geva/peak_ages_V2", full.names = T), value = T, invert = F)

cat.geva <- do.call("rbind",lapply(file.path.geva,
                                        FUN=function(files){
                                          x <- fread(files, header = T)
                                          x$peak <- gsub("data/geva/peak_ages_V2/peak_", "", gsub(".sites.txt","", files))
                                          x
                                        }))

##get random sites
file.path.ran.geva <- grep("*sites.txt",list.files("data/geva/random/array_output", full.names = T), value = T, invert = F)

cat.ran.geva <- do.call("rbind",lapply(file.path.ran.geva,
                                   FUN=function(files){
                                     x <- fread(files, header = T)
                                     x$peak <- "Random"
                                     x
                                   }))

#some exploratory plots
cat.geva.filt <- cat.geva %>% filter(Filtered == 1 & Clock == "J")
cat.ran.geva.filt <- cat.ran.geva %>% filter(Filtered == 1 & Clock == "J")
cat.geva.M <- cat.geva %>% filter(Filtered == 1 & Clock == "M")
cat.ran.geva.M <- cat.ran.geva %>% filter(Filtered == 1 & Clock == "M")

#multiply mean age by generation time
cat.geva.filt$Age <- cat.geva.filt$PostMean*5
cat.ran.geva.filt$Age <- cat.ran.geva.filt$PostMean * 5
cat.geva.M$Age <- cat.geva.M$PostMean*5
cat.ran.geva.M$Age <- cat.ran.geva.M$PostMean*5

cat.geva.sum <- cat.geva.filt %>% 
  group_by(peak) %>% 
  summarise(age = mean(PostMean) * 5)

peak.stats <- read.csv("assem_finch_code/supplementary_data_for_github/Supplementary_data_2_all_peaks_median_diversity_stat.csv")
peak.stats$peak <- as.character(peak.stats$peak)

peak.stats <- left_join(peak.stats, cat.geva.sum, by = "peak")

peak.stats.f <- peak.stats %>% 
  filter(contrast == "fuliginosa_magnirostris") 
peak.stats.f %>% 
  ggplot() + geom_point(aes(x = age, y = time)) +
  geom_abline() +
  xlim(2.5e5, 9e5) +
  ylim(2.5e5, 9e5)+
  labs(x = "geva estimate", y = "dxy estimate")

ggsave("output/geva/compare_dxy_geva.pdf", width = 8, height = 6)

cat.geva.filt %>% 
  filter(N_Concordant > 925) %>% 
  filter(peak == 1) %>% 
  ggplot() + geom_point(aes(x = N_Concordant, y = PostMean))


cat.geva.filt %>% 
  ggplot() + 
  geom_point(aes(x = MarkerID, y = PostMean)) +
  facet_grid(~peak) +
  theme(axis.text.x = element_blank()) 



# ridge plot --------------------------------------------------------------

merge.geva.filt <- rbind(cat.geva.filt, cat.ran.geva.filt)

#merge.geva.filt$peak <- factor(as.numeric(merge.geva.filt$peak))

merge.geva.filt <- merge.geva.filt %>% 
  group_by(peak) %>% 
  summarise(med_time = median(Age),
            cl = quantile(Age, 0.975)) %>% 
  arrange(med_time) %>% 
  mutate(peak_order = 1:n()) %>% 
  select(-med_time) %>% 
  right_join(merge.geva.filt, by = "peak")

#brute force the random SNPs to be on the bottom
merge.geva.filt <- merge.geva.filt %>% 
  mutate(peak_order = ifelse(peak == "Random", 0, peak_order),
         label = as.factor(paste(peak, round(cl/1000000, 2), sep = " ")))

merge.geva.filt$peak <- as.factor(merge.geva.filt$peak)

ggplot(merge.geva.filt, aes(x = Age/1000000, y = fct_reorder(label, peak_order), fill = stat(x))) +
  geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = c(0.975),
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 2, point_alpha = 1) +
  theme_ridges() +
  #scale_fill_cyclical(values = c("grey80", "grey90")) +
  scale_fill_viridis_c(name = "Time", option = "C") +
  labs(y = "Peak", x = "Time (millions of years)")

ggsave("output/geva/geva_time_estimates.pdf", height = 11, width = 8.5)


# estimates for paper -----------------------------------------------------

mean(merge.geva.filt$cl)
sd(merge.geva.filt$cl)
