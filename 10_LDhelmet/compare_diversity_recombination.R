library(tidyverse)
library(data.table)



df.pixy <- read.csv("output/pixy/Four_geospiza_processed_Pixy.csv") %>% 
  select(CHROM, BIN_START, BIN_END, starts_with("pi"))
df.recomb <- read.csv("output/LDHelmet/windowed_rho_20kb.csv")

df.pixy.recomb <- inner_join(df.pixy, df.recomb, by = c("CHROM" = "CHROM", "BIN_START" = "win_start"))

df.pixy.recomb.A <- df.pixy.recomb %>% filter(CHROM!="chrZ")
df.pixy.recomb.Z <- df.pixy.recomb %>% filter(CHROM=="chrZ")

cor.test(df.pixy.recomb.A$pi_fortis, df.pixy.recomb.A$rho_per_kb)
cor.test(df.pixy.recomb.A$pi_magnirostris, df.pixy.recomb.A$rho_per_kb)
cor.test(df.pixy.recomb.A$pi_fuliginosa, df.pixy.recomb.A$rho_per_kb)

cor.test(df.pixy.recomb.Z$pi_fortis, df.pixy.recomb.Z$rho_per_kb)
cor.test(df.pixy.recomb.Z$pi_magnirostris, df.pixy.recomb.Z$rho_per_kb)
cor.test(df.pixy.recomb.Z$pi_fuliginosa, df.pixy.recomb.Z$rho_per_kb)


# fst all speceies --------------------------------------------------------


df.pixy.fst <- read.csv("output/pixy/Four_geospiza_processed_Pixy.csv") %>% 
  select(CHROM, BIN_START, BIN_END, starts_with("fst"))
mean(df.pixy.fst$fst_fortis_fuliginosa, na.rm = T)
mean(df.pixy.fst$fst_fortis_magnirostris, na.rm = T)
mean(df.pixy.fst$fst_fuliginosa_magnirostris, na.rm = T)


# parvulus ----------------------------------------------------------------
# pi ----------------------------------------------------------------------

pixy_contrast_pi <- grep("pi.txt",list.files("data/pixy/parvulus/", full.names = T), value = TRUE)

cat_pixy_pi <- do.call("rbind",lapply(pixy_contrast_pi,
                                      FUN=function(files){
                                        x <- fread(files, skip = 1)
                                        names(x) <- c("pop",	"CHROM",	"BIN_START",	"BIN_END",	"pi",	"no_sites",	"count_diffs",	"count_comparisons",	"count_missing")
                                        x$CHROM <- gsub("LGE22", "29", x$CHROM)
                                        x
                                      }))


df.pixy.recomb2 <- inner_join(cat_pixy_pi, df.recomb, by = c("CHROM" = "CHROM", "BIN_START" = "win_start"))
df.pixy.recomb2.A <- df.pixy.recomb2 %>% filter(CHROM!="chrZ")
df.pixy.recomb2.Z <- df.pixy.recomb2 %>% filter(CHROM=="chrZ")

summary(lm(df.pixy.recomb2.A$pi ~ df.pixy.recomb2.A$rho_per_kb, data = df.pixy.recomb.A))$r.squared
summary(lm(df.pixy.recomb2.Z$pi ~ df.pixy.recomb2.Z$rho_per_kb, data = df.pixy.recomb.Z))$r.squared

plot(df.pixy.recomb2.A$rho_per_kb, df.pixy.recomb2.A$pi)

r.A <- cor(df.pixy.recomb2.A$pi, df.pixy.recomb2.A$rho_per_kb)
#get R2
r.A * r.A
r.Z <- cor(df.pixy.recomb2.Z$pi, df.pixy.recomb2.Z$rho_per_kb)
#get R2
r.Z * r.Z

cor.test(df.pixy.recomb2.A$pi, df.pixy.recomb2.A$rho_per_kb)
cor.test(df.pixy.recomb2.Z$pi, df.pixy.recomb2.Z$rho_per_kb)
