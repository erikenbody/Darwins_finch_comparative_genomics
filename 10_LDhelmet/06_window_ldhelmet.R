library(tidyverse)
library(windowscanr)
library(data.table)
library(patchwork)

file.path.rho <- grep("Icon",list.files("data/LDHelmet/results_run2", full.names = T), value = T, invert = T)
df.rho.files <- data.frame(path = file.path.rho,
                              chrom = gsub("data/LDHelmet/results_run2/","",gsub("_parvulus_ldhelmet.txt","",file.path.rho)))
df.rho.files$chrom <- gsub("LGE22", "29", df.rho.files$chrom)
chr_order3 <- c("chr1", "chr1A" , "chr2", "chr3", "chr4","chr4A", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                "chr11","chr12", "chr13", "chr14", "chr15", "chr17", "chr18", "chr19", "chr20", "chr21", 
                "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chrZ",  "null")

df.rho.files$chr_ordered <- factor(df.rho.files$chrom, levels = chr_order3)
df.rho.files <- df.rho.files %>% arrange(chr_ordered)

# load into one datafile --------------------------------------------------
#file.path.rho <- c("data/LDHelmet/results/chr20_parvulus_ldhelmet.txt","data/LDHelmet/results/chr28_parvulus_ldhelmet.txt")
cat.rho.files <- do.call("rbind",lapply(file.path.rho,
                                       FUN=function(files){
                                         x <- fread(files, skip = 5)
                                         names(x) <- c("left_snp","right_snp", "mean", "p0.025", "p0.500", "p0.098")
                                         chr <- gsub("data/LDHelmet/results_run2/","",gsub("_parvulus_ldhelmet.txt","",files))
                                         x$CHROM <- chr
                                         x$CHROM <- gsub("LGE22", "29", x$CHROM)
                                         x
                                       }))
#head(cat.rho.files)

chr.plots <- list()  
pdf("output/LDHelmet/all_chr_rho_100kb.pdf", width = 8, height = 3)

for (chr in df.rho.files$chrom){
#for (chr in unique(cat.rho.files$CHROM)){
  #df.ldhelment <- fread(chr, skip = 5)
  df.ldhelment <- cat.rho.files %>% filter(CHROM == chr)
  df.ldhelment <- df.ldhelment %>% mutate(distance = right_snp - left_snp,
                          sum_rho = distance * mean)
  df.ldhelment <- df.ldhelment %>% filter(distance < 5000)
  df.windows <- winScan(x = df.ldhelment, 
                          position = "left_snp", 
                          values = "sum_rho", 
                          win_size = 100000,
                          win_step = 100000,
                          funs = "sum")
  
  chr.plots[[chr]] <- ggplot() + 
    geom_line(data = df.windows, aes(x = win_start / 1000000, y = sum_rho_sum / 100, color = sum_rho_sum/100), size = 1.25) +
    theme_bw() + 
    scale_colour_gradient2(low="blue", high="red", mid = "orange", midpoint = (max(df.windows$sum_rho_sum/100, na.rm = TRUE) /2 )) +
    xlab("Position (Mb)") + ylab("Recombination frequency (Rho/kb)") +
    theme(legend.position='none') + ggtitle(chr)
  
  print(chr.plots[[chr]])
  
}
dev.off()

df.ldhelment <- cat.rho.files %>% mutate(distance = right_snp - left_snp,
                                        sum_rho = distance * mean)

df.ldhelment <- df.ldhelment %>% filter(distance < 5000)

df.windows <- winScan(x = df.ldhelment, 
                      position = "left_snp", 
                      values = "sum_rho", 
                      win_size = 20000,
                      win_step = 20000,
                      groups = "CHROM",
                      funs = "sum")
df.windows <- df.windows %>% mutate(rho_per_kb = sum_rho_sum/20)
write_csv(df.windows,"output/LDHelmet/windowed_rho_20kb.csv")

df.windows100 <- winScan(x = df.ldhelment, 
                      position = "left_snp", 
                      values = "sum_rho", 
                      win_size = 100000,
                      win_step = 100000,
                      groups = "CHROM",
                      funs = "sum")
df.windows100 <- df.windows100 %>% mutate(rho_per_kb = sum_rho_sum/100)

write_csv(df.windows100,"output/LDHelmet/windowed_rho_100kb.csv")

# centimorgan simplified --------------------------------------------------

cat.rho.files.convert <- cat.rho.files %>% mutate(distance = right_snp - left_snp,
                                                  r = (mean / (4 * 33861)) , #calc cM per Mb by converting with Ne (from Nature 2015 finch paper)
                                                  rcM = (r * 100)) %>%  #using morgan function
  group_by(CHROM) %>% mutate(cM = cumsum(rcM *distance)) #cumulatively sum of these values to get an estimate of cM

chr.tab.cm <- cat.rho.files.convert %>% group_by(CHROM) %>% 
  summarise(max_cm = max(cM, na.rm = T)) %>% as.data.frame()
sum(chr.tab.cm$max_cm)

chr.plots <- list()  

for (chr in df.rho.files$chrom){
  #df.ldhelment <- fread(chr, skip = 5)
  df.ldhelment <- cat.rho.files.convert %>% filter(CHROM == chr)
  
  df.windows <- winScan(x = df.ldhelment, 
                        position = "left_snp", 
                        values = "rcM", 
                        win_size = 100000,
                        win_step = 100000,
                        funs = "mean")
  
  df.windows$cMMb.roll <- zoo::rollmean(df.windows$rcM_mean, 10, fill = NA, align = "right")
  
  chr.plots[[chr]] <- ggplot() + 
    geom_line(data = df.windows, aes(x = win_start / 1000000, y = cMMb.roll*1000000), size = 1.25) +
    theme_bw() + 
    #scale_colour_gradient2(low="blue", high="red", mid = "orange", midpoint = (max(df.windows$sum_rho_sum/100, na.rm = TRUE) /2 )) +
    xlab("Position (Mb)") + ylab("cM / Mb") +
    theme(legend.position='none') + ggtitle(chr)
  
}

wrap_plots(chr.plots,ncol = 6)
ggsave("output/LDHelmet/wrapped_cMMb_plots_V2.png", height = 10, width = 12)


#output a centimorgan file

for (chr in unique(cat.rho.files.convert$CHROM)){
  df.sub.out <- cat.rho.files.convert %>% 
    ungroup() %>% 
    filter(CHROM == chr) %>% 
    mutate(cMMb = rcM * 1000000) %>% 
    select(left_snp, cMMb, cM)
  names(df.sub.out) <- c("Position(bp)", "Rate(cM/Mb)",	"Map(cM)")
  write_tsv(df.sub.out, paste0("output/LDHelmet/recombination_rate_per_chr_LDHELMET/",chr,"_recomb_rate.txt"), col_names = T)
}


