library(data.table)
library(tidyverse)

source("assem_finch_code/07_popgen/generic_manhattan_fnc.R")

fortis.files <- grep("fortis",list.files("data/tajimas_d/with_invariant", full.names = T, recursive = T), value = TRUE)
fortis.taj <- do.call("rbind",lapply(fortis.files,
                                      FUN=function(files){
                                        x <- fread(files, skip = 1)
                                        names(x) <- c("CHROM",	"BIN_START", "N_SNPS", "tajd_fortis")
                                        x
                                      }))
#

fuliginosa.files <- grep("fuliginosa",list.files("data/tajimas_d/with_invariant", full.names = T, recursive = T), value = TRUE)
fuliginosa.taj <- do.call("rbind",lapply(fuliginosa.files,
                                     FUN=function(files){
                                       x <- fread(files, skip = 1)
                                       names(x) <- c("CHROM",	"BIN_START", "N_SNPS", "tajd_fuliginosa")
                                       x
                                     }))

#


magnirostris.files <- grep("magnirostris",list.files("data/tajimas_d/with_invariant", full.names = T, recursive = T), value = TRUE)
magnirostris.taj <- do.call("rbind",lapply(magnirostris.files,
                                     FUN=function(files){
                                       x <- fread(files, skip = 1)
                                       names(x) <- c("CHROM",	"BIN_START", "N_SNPS", "tajd_magnirostris")
                                       x
                                     }))


# -------------------------------------------------------------------------


df.taj.d <- cbind(fortis.taj, fuliginosa.taj[,"tajd_fuliginosa"], magnirostris.taj[,"tajd_magnirostris"]) %>% 
  mutate(BIN_END = BIN_START + 19999,
         BIN_MID = BIN_START + ((BIN_END - BIN_START) /2))
df.taj.d$CHROM <- factor(df.taj.d$CHROM, levels = chr_order)
df.taj.d <- df.taj.d %>% arrange(CHROM)
write_csv(df.taj.d, "data/tajimas_d/with_invariant/three_species_tajimasd_with_invariant.csv")
