suppressMessages(library(zoo))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(R.utils))
suppressMessages(library(ggrepel))
suppressMessages(library(patchwork))


# -------------------------------------------------------------------------

chr_order <- c("chr1", "chr1A" , "chr2", "chr3", "chr4","chr4A", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
               "chr11","chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
               "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chrZ",  "chrunknown")

# FUNCTIONS ---------------------------------------------------------------
#fst.path <- "data/fst_scans/SGF_MGF.windowed.weir.fst"
fst.order <- function(fst.path){
  fst <- fread(fst.path)
  fst$chr_ordered <- factor(fst$CHROM, levels = chr_order)
  fst <- fst %>% dplyr::arrange(chr_ordered, BIN_START)
  fst$row<-1:nrow(fst)
  fst$fstrollmean <- zoo::rollmean(fst$WEIGHTED_FST,50,fill=NA)
  fst$chr_labels <- gsub("chr", "", fst$chr_ordered)
  chr_breaks <- fst %>% group_by(chr_labels) %>% dplyr::summarise(chr_breaks = mean(row))
  fst
  

}

zfst <- function(df.in){
  zdf <- df.in %>% filter(CHROM == "chrZ") %>% mutate(zfst = (WEIGHTED_FST - mean(WEIGHTED_FST, na.rm = T)) / sd(WEIGHTED_FST, na.rm =T))
  adf <- df.in %>% filter(CHROM != "chrZ" | is.na(CHROM)) %>% mutate(zfst = (WEIGHTED_FST - mean(WEIGHTED_FST, na.rm = T)) / sd(WEIGHTED_FST, na.rm = T))
  fst.comb <- rbind(adf, zdf)
  fst.comb
}

zfst2 <- function(df.in, input.var){
  zdf <- df.in %>% filter(CHROM == "chrZ") %>% mutate(zfst = (input.var - mean(input.var, na.rm = T)) / sd(input.var, na.rm =T))
  adf <- df.in %>% filter(CHROM != "chrZ" | is.na(CHROM)) %>% mutate(zfst = (input.var - mean(input.var, na.rm = T)) / sd(input.var, na.rm = T))
  fst.comb <- rbind(adf, zdf)
  fst.comb
}

#df.in <- fst.order(fst.path)
#zfst(df.in)

#df.in <- df.comp.badmap
#typical manhattan
manc <- function(df.in, input.var){
  chr_breaks <- df.in %>% filter(chr_ordered != "chrunknown" & !is.na(row)) %>% 
    mutate(chr_ordered = factor(chr_ordered, levels = chr_order)) %>%
    group_by(chr_ordered, chr_labels) %>% 
    dplyr::summarise(chr_breaks = mean(row), .groups = "keep")
  
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
          #axis.text.x = element_text(angle = 45, color = "black"),
          #axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.y = element_text(size=10),
          axis.text = element_text(size=10),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    geom_point(size=0.9,shape=20,stroke=0.2) +
    scale_color_manual(values=rep(c("grey30","grey70"))) +
    #scale_color_manual(values = c(rep_len(c("grey30", "red"), length(unique(chr_breaks$chr_ordered))+1))) #+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                       labels = function(labels) {
                         sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                       }) +
    labs(y=expression(F[ST])) 
    #scale_x_continuous(breaks=chr_breaks$chr_breaks, 
    #                   labels = chr_breaks$chr_labels)
}

#peak_list <- function(input.df, input.var){
#  #input.df$var <- as.list(input.df[,input.var])
#  input.df <- as.data.frame(input.df)
#  #fst <- fst %>% filter(chr_ordered == "chr1A")
#  auto_fst <- input.df %>% filter(chr_ordered != "chrZ" & chr_ordered!="chrunknown")
#  Z_fst <- input.df %>% filter(chr_ordered == "chrZ")
#  
#  #background.chr <- input.df %>% filter(chr_ordered == "chr10" & chr.pos < 15000)
#  
#  #ggplot(background.chr) + geom_point(aes(x = chr.pos, y = LR))
#  
#  #peaks_fst_auto <- auto_fst %>% filter(fstrollmean > quantile(auto_fst$fstrollmean,0.95,na.rm=T)[[1]]) 
#  peaks_fst_auto <- auto_fst %>% filter(!!sym(input.var) > quantile(auto_fst[,input.var],0.995,na.rm=T)[[1]]) #was 0.989
#  #peaks_fst_auto <- auto_fst %>% filter(!!sym(input.var) > quantile(auto_fst[,input.var],0.995,na.rm=T)[[1]]) 
#  
#  peaks_fst_Z <- Z_fst %>% filter(!!sym(input.var) > quantile(Z_fst[,input.var],0.995,na.rm=T)[[1]]) 
#  peaks_fst <- rbind(peaks_fst_auto, peaks_fst_Z)
#  
#  fst_range <- GRanges(peaks_fst$CHROM, IRanges(as.numeric(peaks_fst$BIN_START), as.numeric(peaks_fst$BIN_END)))
#  
#  range.peaks_fst <- GRanges(peaks_fst$CHROM, IRanges(as.numeric(peaks_fst$BIN_START), as.numeric(peaks_fst$BIN_END)))
#  reduce.range.peaks_fst <- GenomicRanges::reduce(range.peaks_fst, min.gapwidth = 1000) #merge closest window
#  names(reduce.range.peaks_fst) <- 1:length(reduce.range.peaks_fst)
#  
#  datalist <- list()
#  for (i in (1:length(reduce.range.peaks_fst))){
#    x<-peaks_fst[queryHits(findOverlaps(fst_range, reduce.range.peaks_fst[i, ])), ]
#    x$peak <- i #why the error here?
#    datalist[[i]] <- x # add it to your list
#  }
#  peaks.named = do.call(rbind, datalist) 
#  
#  f.peaks.named.A <- peaks.named %>% filter(chr_ordered!="chrZ") %>% 
#    group_by(peak) %>% 
#    filter(max(!!sym(input.var)) > quantile(auto_fst[,input.var],0.995,na.rm=T)[[1]]) 
#  
#  f.peaks.named.Z <- peaks.named %>% filter(chr_ordered=="chrZ") %>% 
#    group_by(peak) %>% 
#    filter(max(!!sym(input.var)) > quantile(Z_fst[,input.var],0.995,na.rm=T)[[1]]) 
#  
#  f.peaks.named <- rbind(f.peaks.named.A, f.peaks.named.Z)
#  
#  f.peaks.named$peak2 <- f.peaks.named %>% ungroup() %>% group_indices(peak)
#  f.peaks.named$peak <- f.peaks.named$peak2 ; f.peaks.named$peak2 <- NULL
#  f.peaks.named
#  
#}

peak_list <- function(input.df, input.var, threshold_1, threshold_2){
  #input.df$var <- as.list(input.df[,input.var])
  input.df <- as.data.frame(input.df)
  
  #fst <- fst %>% filter(chr_ordered == "chr1A")
  auto_fst <- input.df %>% filter(chr_ordered != "chrZ" & chr_ordered!="chrunknown")
  Z_fst <- input.df %>% filter(chr_ordered == "chrZ")
  
  #background.chr <- input.df %>% filter(chr_ordered == "chr10" & chr.pos < 15000)
  
  #ggplot(background.chr) + geom_point(aes(x = chr.pos, y = LR))
  
  #peaks_fst_auto <- auto_fst %>% filter(fstrollmean > quantile(auto_fst$fstrollmean,0.95,na.rm=T)[[1]]) 
  peaks_fst_auto <- auto_fst %>% filter(!!sym(input.var) > quantile(auto_fst[,input.var],threshold_1,na.rm=T)[[1]]) 
  #peaks_fst_auto <- auto_fst %>% filter(!!sym(input.var) > quantile(auto_fst[,input.var],0.995,na.rm=T)[[1]]) 
  
  peaks_fst_Z <- Z_fst %>% filter(!!sym(input.var) > quantile(Z_fst[,input.var],threshold_1,na.rm=T)[[1]]) 
  peaks_fst <- rbind(peaks_fst_auto, peaks_fst_Z)
  
  fst_range <- GRanges(peaks_fst$CHROM, IRanges(as.numeric(peaks_fst$BIN_START), as.numeric(peaks_fst$BIN_END)))
  
  range.peaks_fst <- GRanges(peaks_fst$CHROM, IRanges(as.numeric(peaks_fst$BIN_START), as.numeric(peaks_fst$BIN_END)))
  reduce.range.peaks_fst <- GenomicRanges::reduce(range.peaks_fst, min.gapwidth = 100000) #merge closest window
  names(reduce.range.peaks_fst) <- 1:length(reduce.range.peaks_fst)
  
  datalist <- list()
  for (i in (1:length(reduce.range.peaks_fst))){
    x<-peaks_fst[queryHits(findOverlaps(fst_range, reduce.range.peaks_fst[i, ])), ]
    x$peak <- i #why the error here?
    datalist[[i]] <- x # add it to your list
  }
  peaks.named = do.call(rbind, datalist) 
  
  f.peaks.named.A <- peaks.named %>% filter(chr_ordered!="chrZ") %>% 
    group_by(peak) %>% 
    filter(max(!!sym(input.var)) > quantile(auto_fst[,input.var],threshold_2,na.rm=T)[[1]]) 
  
  f.peaks.named.Z <- peaks.named %>% filter(chr_ordered=="chrZ") %>% 
    group_by(peak) %>% 
    filter(max(!!sym(input.var)) > quantile(Z_fst[,input.var],threshold_2,na.rm=T)[[1]]) 
  
  f.peaks.named <- rbind(f.peaks.named.A, f.peaks.named.Z)
  
  f.peaks.named$peak2 <- f.peaks.named %>% ungroup() %>% group_indices(peak)
  f.peaks.named$peak <- f.peaks.named$peak2 ; f.peaks.named$peak2 <- NULL
  f.peaks.named
  
}

peak_bed <- function(peaks){
  peak.bed <- peaks %>% group_by(peak, CHROM) %>% 
    summarise(start = min(BIN_START),
              end = max(BIN_END))
}

peak_bed_row <- function(peaks){
  peak.bed <- peaks %>% group_by(peak, CHROM) %>% 
    summarise(start = min(row),
              end = max(row))
}

#peaks with diversity also marked
plot_peaks2 <- function(comp.order, peak_file, input.var){
  comp.order <- as.data.frame(comp.order)
  chr_breaks <- comp.order %>% filter(chr_ordered != "chrunknown" & !is.na(row)) %>% 
    mutate(chr_ordered = factor(chr_ordered, levels = chr_order)) %>%
    group_by(chr_ordered, chr_labels) %>% 
    dplyr::summarise(chr_breaks = mean(row))
  
  chrom.colors <- data.frame(chr_ordered = grep("chr", unique(comp.order$chr_ordered), value = T),
                             color.num = rep(1:2,length(grep("chr", unique(comp.order$chr_ordered))))) %>% 
    distinct(chr_ordered, .keep_all = T)
  
  df.in2 <- comp.order %>% #mutate(row = 1:n()) %>% 
    left_join(chrom.colors, by = "chr_ordered") %>% 
    mutate(color.num = as.factor(color.num))
  
  max.val <- max(comp.order[,input.var], na.rm = T)
  
  df.in2 %>% filter(chr_labels != "unknown" & !is.na(row)) %>%
    ggplot(aes_string(x = "row", y = input.var, col = "color.num")) + theme_bw() +
    theme(legend.position="none",
          panel.border=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, color = "black"),
          panel.grid = element_blank(),
          panel.grid.major.y=element_line(color="grey60",size=0.2),
          panel.grid.minor.y=element_line(color="grey60",size=0.1),
          axis.title.y = element_text(size=20),
          axis.text = element_text(size=12),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    #geom_point(data = subset(df.in2, diversity == "low"), aes_string(x = "row", y = input.var, col = "color.num"), size = 1.5, color = "gold1")+
    geom_point(size=0.9,shape=20,stroke=0.2) +
    scale_color_manual(values=rep(c("grey30","grey58"))) +
    #scale_color_manual(values = c(rep_len(c("grey30", "red"), length(unique(chr_breaks$chr_ordered))+1))) #+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    labs(y=expression(F[ST])) +
    #scale_x_continuous(breaks=chr_breaks$chr_breaks, 
    #                   labels = chr_breaks$chr_labels) +
    #https://stackoverflow.com/questions/50399838/how-to-alternate-a-new-line-for-overlapping-x-axis-labels
    scale_x_continuous(breaks = chr_breaks$chr_breaks, 
                       labels = function(labels) {
                         sapply(seq_along(labels), function(i) paste0(ifelse(i %% 2 == 0, '', '\n'), chr_breaks$chr_labels[i]))
                       }) +
    geom_point(data = peak_file, aes(x = start, y = max.val + 0.02), color = "blue", alpha = 0.4) 
  
}

#comp.order <- df.sgf.mgf
#comp_row.bed <- df.sgf.mgf.peaks.bed.row
#df.sgf.mgf.peaks.genes -> comp.annotated
#chr.sub <- "chr4"
#input.var <- "zfst"

plot_peaks_subset <-function(comp.order, comp_row.bed, comp.annotated, chr.sub, input.var){

  
  df.subset <- subset(comp.order, chr_ordered == chr.sub)
  
  df.dummy <- comp.order %>% select(chr_ordered, CHROM) %>% 
    distinct(CHROM, .keep_all = T) %>% 
    mutate(start = NA, start.y = NA, end = NA, gene_name = NA, peak = 9999) %>% 
    filter(!CHROM %in% comp_row.bed$CHROM) 
  
  #comp_row.bed <- bind_rows(comp_row.bed, df.dummy[c(2,3,5,7)])
  comp_row.bed <- bind_rows(comp_row.bed, df.dummy)
  #ann.label <- ann.label %>% full_join(comp_row.bed, by = c("CHROM"))
  no_genes <- comp_row.bed %>% filter(!CHROM %in% comp.annotated$CHROM) %>% mutate(gene_name = "(no genes)")
  ann.label <- bind_rows(comp.annotated, no_genes)
  head(comp_row.bed)
  head(ann.label)
  
  #df.x <- comp_row.bed %>% select(peak, start) %>% dplyr::rename(start.peak = start)
  ann.label <- ann.label %>% filter(CHROM %in% df.subset$CHROM) %>% select(-start.y) %>%  left_join(comp_row.bed[,-c(5,6,7)], by = c("peak", "CHROM"))
  df.subset<-as.data.frame(df.subset)
  max.val <- max(df.subset[,input.var], na.rm = T)
  
  #comp.order$chr.pos.bp <- comp.order$chr.pos * 5000 / 1000000
  
  p1 <- ggplot(data = base::subset(comp.order, chr_ordered == chr.sub), aes_string(x="row",y=input.var)) + theme_bw() +
    theme(legend.position="none",
          #panel.border=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          #panel.grid = element_blank(),
          panel.grid.major.y=element_line(color="grey60",size=0.2),
          panel.grid.minor.y=element_line(color="grey60",size=0.1),
          axis.title.y = element_text(size=20),
          axis.text = element_text(size=12),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    #geom_point(data = subset(comp.order, diversity == "low" & chr_ordered == chr.sub), aes_string(x="row",y=input.var), size = 3, color = "gold1")+
    geom_point(size=1.1,shape=19) +
    labs(y = expression(ZF[ST]), title = chr.sub) +
    geom_segment(data = base::subset(comp_row.bed, CHROM %in% df.subset$CHROM), aes(y = 16, yend = 16, x = start - 50, xend = end + 50), color = "blue", size = 2) +
    geom_label_repel(data = base::subset(ann.label, CHROM %in% df.subset$CHROM), aes(x = start.y, y = 16, label = gene_name))
  
  p1 
  
}

## March 3 update
plot_peaks_subset2 <-function(comp.order, comp.bed, comp.annotated, chr.sub, input.var){
  
  
  df.subset <- subset(comp.order, chr_ordered == chr.sub)
  
  df.dummy <- comp.order %>% select(chr_ordered, CHROM) %>% 
    distinct(CHROM, .keep_all = T) %>% 
    mutate(start = NA, start.y = NA, end = NA, gene_name = NA, peak = 9999) %>% 
    filter(!CHROM %in% comp.bed$CHROM) 
  
  #comp.bed <- bind_rows(comp.bed, df.dummy[c(2,3,5,7)])
  comp.bed <- bind_rows(comp.bed, df.dummy)
  #ann.label <- ann.label %>% full_join(comp.bed, by = c("CHROM"))
  no_genes <- comp.bed %>% filter(!CHROM %in% comp.annotated$CHROM) %>% mutate(gene_name = "(no genes)")
  ann.label <- bind_rows(comp.annotated, no_genes)
  head(comp.bed)
  head(ann.label)
  
  #df.x <- comp.bed %>% select(peak, start) %>% dplyr::rename(start.peak = start)
  ann.label <- ann.label %>% filter(CHROM %in% df.subset$CHROM) %>% select(-start.y) %>%left_join(comp.bed[,-c(5:7)], by = c("peak", "CHROM"))
  
  df.subset<-as.data.frame(df.subset)
  max.val <- max(df.subset[,input.var], na.rm = T)
  
  #comp.order$chr.pos.bp <- comp.order$chr.pos * 5000 / 1000000
  
  p1 <- ggplot(data = base::subset(comp.order, chr_ordered == chr.sub), aes_string(x="BIN_START",y=input.var)) + theme_bw() +
    theme(legend.position="none",
          panel.border=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major.y=element_line(color="grey60",size=0.2),
          panel.grid.minor.y=element_line(color="grey60",size=0.1),
          axis.title.y = element_text(size=20),
          axis.text = element_text(size=12),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_line(size=0.2)) +
    #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1),minor_breaks = NULL) +
    #geom_point(data = subset(comp.order, diversity == "low" & chr_ordered == chr.sub), aes_string(x="row",y=input.var), size = 3, color = "gold1")+
    geom_point(size=1.1,shape=19) +
    labs(y = expression(ZF[ST]), title = chr.sub) +
    #labs(y = input.var, title = chr.sub) +
    geom_segment(data = base::subset(comp.bed, CHROM %in% df.subset$CHROM), aes(y = 10, yend = 10, x = start - 50, xend = end + 50), color = "blue", size = 2) +
    geom_label_repel(data = base::subset(ann.label, CHROM %in% df.subset$CHROM), aes(x = start.y, y = 11, label = gene_name), max.overlaps = 50) +
    scale_x_continuous(labels = function(l) {
      paste0(round(l/1e6,1),"Mbp")
    })  
  
  p1 
  
}


input_pixy_V2 <- function(pixy_contrast){
  #pixy_contrast <- "data/pixy/subspecies/Mot_fla_fla_Mot_fla_thu/"
  pixy_contrast <- pixy_contrast
  # pi ----------------------------------------------------------------------
  
  pixy_contrast_pi <- grep("pi.txt",list.files(pixy_contrast, full.names = T), value = TRUE)
  
  cat_pixy_pi <- do.call("rbind",lapply(pixy_contrast_pi,
                                        FUN=function(files){
                                          x <- fread(files, skip = 1)
                                          names(x) <- c("pop",	"CHROM",	"BIN_START",	"BIN_END",	"pi",	"no_sites",	"count_diffs",	"count_comparisons",	"count_missing")
                                          x$CHROM <- gsub("LGE22", "29", x$CHROM)
                                          x
                                        }))
  #extract no_sites, which is same in all populations, to be added to the main dataframe
  
  no_sites.df <- cat_pixy_pi %>% filter(pop == unique(cat_pixy_pi$pop)[1]) %>% select(CHROM, BIN_START, BIN_END, no_sites)
  
  cat_pixy_pi.name <- cat_pixy_pi %>% 
    mutate(pop_num = paste0("pi_", pop))
  
  #contrast_name <- paste0(unique(cat_pixy_pi.name$pop), collapse = "_")
  
  cat_pixy_pi.long <- cat_pixy_pi.name %>%   
    select(-pop, -no_sites, -count_diffs, -count_comparisons, -count_missing)  %>% 
    pivot_wider(
      names_from = pop_num,
      names_sep = "_",
      values_from = pi
    )
  
  cat_pixy_pi.long <- left_join(cat_pixy_pi.long, no_sites.df, by = c("CHROM","BIN_START","BIN_END"))
  
  rm(cat_pixy_pi.name)
  
  # dxy ---------------------------------------------------------------------
  
  pixy_contrast_dxy <- grep("dxy.txt",list.files(pixy_contrast, full.names = T), value = TRUE)
  
  cat_pixy_dxy <- do.call("rbind",lapply(pixy_contrast_dxy,
                                         FUN=function(files){
                                           x <- fread(files, skip = 1)
                                           names(x) <- c("pop1",	"pop2",	"CHROM", "BIN_START",	"BIN_END",	"dxy", "no_sites", "count_diffs",	"count_comparisons",	"count_missing")
                                           x$CHROM <- gsub("LGE22", "29", x$CHROM)
                                           x
                                         }))
  
  cat_pixy_dxy <- cat_pixy_dxy %>% mutate(contrast = paste0("dxy_", pop1, "_", pop2))
  
  cat_pixy_dxy.long <- cat_pixy_dxy %>%   
    select(-pop1, -pop2, -no_sites, -count_diffs, -count_comparisons, -count_missing)  %>% 
    pivot_wider(
      names_from = contrast,
      names_sep = "_",
      values_from = dxy
    )
  
  rm(cat_pixy_dxy)
  
  #same number of sites in all comparisons
  #library(pysch)
  #no.sites <- cat_pixy_dxy.long %>% select(contains("no_sites"))
  #corr.test(no.sites)
  
  
  # fst ---------------------------------------------------------------------
  
  pixy_contrast_fst <- grep("fst.txt",list.files(pixy_contrast, full.names = T), value = TRUE)
  
  cat_pixy_fst <- do.call("rbind",lapply(pixy_contrast_fst,
                                         FUN=function(files){
                                           x <- fread(files, skip = 1)
                                           names(x) <- c("pop1",	"pop2",	"CHROM", "BIN_START",	"BIN_END",	"WEIGHTED_FST", "SNPs")
                                           x$CHROM <- gsub("LGE22", "29", x$CHROM)
                                           x
                                         }))
  
  #reset fst to 0 if its <0
  cat_pixy_fst<-cat_pixy_fst %>% mutate(WEIGHTED_FST = ifelse(WEIGHTED_FST < 0, 0, WEIGHTED_FST))
  
  cat_pixy_fst <- cat_pixy_fst %>% mutate(contrast = paste0("fst_", pop1, "_", pop2))
  
  cat_pixy_fst.long <- cat_pixy_fst %>%   
    select(-pop1, -pop2)  %>% 
    pivot_wider(
      names_from = contrast,
      names_sep = "_",
      values_from = WEIGHTED_FST
    )
  rm(cat_pixy_fst)
  
  # merge -------------------------------------------------------------------
  
  df.merge1 <- full_join(cat_pixy_dxy.long, cat_pixy_fst.long, by = c("CHROM","BIN_START","BIN_END"))
  df.merge2 <- full_join(df.merge1, cat_pixy_pi.long, by = c("CHROM","BIN_START","BIN_END"))
  
  #no.sites <- df.merge2 %>% select(contains("no_sites"))
  #corr.test(no.sites)
  
  df.pixy_out <- df.merge2# %>% select(CHROM, BIN_START, BIN_END, dxy, WEIGHTED_FST, contains("pi"), no_sites, SNPs)
  
  df.pixy_out$chr_ordered <- factor(df.pixy_out$CHROM, levels = chr_order)
  df.pixy_out <- df.pixy_out %>% dplyr::arrange(chr_ordered, BIN_START)
  
  df.pixy_out$row<-1:nrow(df.pixy_out)
  #df.pixy_out$fstrollmean <- zoo::rollmean(df.pixy_out$WEIGHTED_FST,50,fill=NA)
  #df.pixy_out$dxyrollmean <- zoo::rollmean(df.pixy_out$dxy,50,fill=NA)
  #df.pixy_out$pi_pop1rollmean <- zoo::rollmean(df.pixy_out$pi_pop1,50,fill=NA)
  #df.pixy_out$pi_pop2rollmean <- zoo::rollmean(df.pixy_out$pi_pop2,50,fill=NA)
  #df.pixy_out<-zfst(df.pixy_out)
  
  df.pixy_out$chr_labels <- gsub("chr", "", df.pixy_out$chr_ordered)
  
  
  
  return(df.pixy_out)
}  

plot_pixy_V2 <- function(df.pixy_out, contrast_name){
  library(stringr)
  
  
  contr.spl <- str_split(contrast_name, "_", simplify = T)
  
  pop1 <- paste(contr.spl[1:1], collapse="_") 
  pop2 <- paste(contr.spl[2:2], collapse="_") 
  
  df.pixy_contrast <- df.pixy_out %>% 
    select(CHROM, chr_ordered,chr_labels, BIN_START, BIN_END, row, no_sites,ends_with(contrast_name), paste0("pi_",pop1), paste0("pi_",pop2))
  
  df.pixy_contrast[,"WEIGHTED_FST"] <- df.pixy_contrast[,paste0("fst_",contrast_name)] ;df.pixy_contrast[,paste0("fst_",contrast_name)] <- NULL 
  df.pixy_contrast[,"dxy"] <- df.pixy_contrast[,paste0("dxy_",contrast_name)] ;df.pixy_contrast[,paste0("dxy_",contrast_name)] <- NULL 
  df.pixy_contrast<-zfst(df.pixy_contrast)
  
  # save intermediate output ------------------------------------------------
  
  comp.peak <- peak_list(df.pixy_contrast, "zfst",0.995,0.999)
  fwrite(comp.peak, paste0("output/pixy/",contrast_name, "_peaks.txt"), sep = "\t", quote = F)
  comp.bed <- peak_bed(comp.peak)
  fwrite(comp.bed, paste0("output/pixy/",contrast_name, "_peaks.bed"), sep = "\t", quote = F)
  comp_row.bed <- peak_bed_row(comp.peak)
  fwrite(comp_row.bed, paste0("output/pixy/",contrast_name, "_row_peaks.bed"), sep = "\t", quote = F)
  comp.annotated <- get_genes(comp.bed) 
  fwrite(comp.annotated, paste0("output/pixy/",contrast_name, "_genes.txt"), sep = "\t", quote = F)
  
  
  chr.plots.pixy2 <- list()
  pdf(paste0("output/pixy/",contrast_name,"_pixy_fst_dxy_pi_plot_wide.pdf"), width = 8.5, height = 11, useDingbats = F)
  
  for (chr.sub in grep("chr",unique(df.pixy_contrast$CHROM), value = T)){
    print(chr.sub)
    df.pixy.sub <- df.pixy_contrast %>% filter(CHROM == chr.sub)
    #df.box <- df.low.div.ends.bed %>% filter(CHROM == chr.sub)
    xlim.set <- c(min(df.pixy.sub$BIN_START/1000000), max(df.pixy.sub$BIN_START/1000000))
    #p.zfst <- plot_peaks_subset(df.pixy_contrast, df.sgf.mgf.peaks.bed.row, df.sgf.mgf.peaks.genes, chr.sub, "zfst")
    ymax.set.pi <- max(df.pixy.sub[,paste0("pi_",pop1)], na.rm = T)
    ymax.set.dxy <- max(df.pixy.sub[,"dxy"], na.rm = T)
    
    
    #df.pixy.sub %>% ggplot() + geom_point(aes(x = BIN_START/1000000, y = zfst)) + theme_bw()
    
    #get a little lazy with it
    
    df.pixy.sub$fstrollmean <- zoo::rollmean(df.pixy.sub$WEIGHTED_FST,50,fill=NA)
    df.pixy.sub$dxyrollmean <- zoo::rollmean(df.pixy.sub$dxy,50,fill=NA)
    
    df.pixy.sub$pi_pop1rollmean <- zoo::rollmean(df.pixy.sub[,paste0("pi_",pop1)],50, fill=NA)
    df.pixy.sub$pi_pop2rollmean <- zoo::rollmean(df.pixy.sub[,paste0("pi_",pop2)],50, fill=NA)
    
    ymax.set.zfst <- max(df.pixy.sub$zfst, na.rm = T)
    
    p.zfst <-  plot_peaks_subset2(df.pixy.sub, comp.bed, comp.annotated, chr.sub, "zfst")
    
    p.fst <- df.pixy.sub %>% ggplot() + 
      geom_point(aes(x = BIN_START/1000000, y = WEIGHTED_FST)) + 
      geom_line(aes(x = BIN_START/1000000, y = fstrollmean), color = "red") +
      theme_bw() + labs(y = expression(F[ST]), x = "Mbp") +
      theme(legend.position="none",
            #panel.border=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x = element_blank(),
            #panel.grid = element_blank(),
            panel.grid.major.y=element_line(color="grey60",size=0.2),
            panel.grid.minor.y=element_line(color="grey60",size=0.1),
            axis.title.y = element_text(size=20),
            axis.text = element_text(size=12),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_line(size=0.2)) +
      #geom_rect(data = df.box, mapping = aes(xmin = start/1000000, xmax = end/1000000, ymin = 0, ymax = 1, fill = type), alpha = 0.5) +
      ylim(0,1) + ggtitle(chr.sub)
    
    p.dxy <- df.pixy.sub %>% ggplot() + 
      geom_point(aes(x = BIN_START/1000000, y = dxy)) + 
      geom_line(aes(x = BIN_START/1000000, y = dxyrollmean), color = "red") +
      theme_bw() + labs(y = expression(D[XY]), x = "Mbp") +
      theme(legend.position="none",
            #panel.border=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x = element_blank(),
            #panel.grid = element_blank(),
            panel.grid.major.y=element_line(color="grey60",size=0.2),
            panel.grid.minor.y=element_line(color="grey60",size=0.1),
            axis.title.y = element_text(size=20),
            axis.text = element_text(size=12),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_line(size=0.2)) +
      ylim(0,ymax.set.dxy)
    
    p.pi <- df.pixy.sub %>% ggplot() + 
      #geom_point(aes(x = BIN_START/1000000, y = dxy)) + 
      geom_line(aes(x = BIN_START/1000000, y = pi_pop1rollmean), color = "red", size = 1) +
      geom_line(aes(x = BIN_START/1000000, y = pi_pop2rollmean), color = "blue", size = 1) +
      theme_bw() + labs(y = expression(pi), x = "Mbp") +
      theme(legend.position="none",
            #panel.border=element_blank(),
            axis.title.x=element_blank(),
            #axis.text.x = element_blank(),
            #panel.grid = element_blank(),
            panel.grid.major.y=element_line(color="grey60",size=0.2),
            panel.grid.minor.y=element_line(color="grey60",size=0.1),
            axis.title.y = element_text(size=20),
            axis.text = element_text(size=12),
            #axis.ticks.x=element_blank(),
            axis.ticks.y=element_line(size=0.2)) +
      ylim(0,ymax.set.pi)
    
    #chr.plots.pixy2[[chr.sub]] <- p.zfst / p.fst / p.dxy / p.pi
    chr.plots.pixy2[[chr.sub]] <- p.zfst / p.fst / p.dxy / p.pi
    print(chr.plots.pixy2[[chr.sub]])
  }
  dev.off()
}


# annotation ---------------------------------------------------
#library(rtracklayer)

#wag3.gff <- readGFF("data/annotation/Camarhynchus_parvulus.Camarhynchus_parvulus_V1.1.102.gtf")
#wag3.df <- as.data.frame(wag3.gff)
#
#wag3.df.an <- wag3.df %>% 
#  filter(type == "gene") %>% 
#  select(seqid, start, end, gene_id)
#
#parv.annie <- fread("data/annotation/Camarhynchus_parvulus_V1.1.annie", header = F) %>% 
#  filter(V2 == "name") %>% 
#  mutate(V1 = gsub("gene:","",V1, fixed = T)) %>% 
#  select(-V2, gene_id = V1, Name = V3) %>% 
#  separate(Name, into = c("Name", NULL), sep ="_") %>% 
#  distinct(gene_id, .keep_all = T)
#
#parv.df.an <- left_join(wag3.df.an, parv.annie, by = "gene_id") %>% 
#  mutate(gene_symbol = ifelse(is.na(Name), gene_id, Name))
#
#
#
#write.csv(parv.df.an, "data/annotation/Camarhynchus_parvulus_V1.1_genes.csv", row.names = F)
  
#parv.df.an <- read.csv("data/annotation/Camarhynchus_parvulus_V1.1_genes.csv") %>% 
#  mutate(seqid = paste0("chr",seqid))
#write.csv(parv.df.an, "data/annotation/Camarhynchus_parvulus_V1.1_genes.csv", row.names = F)
parv.df.an <- read.csv("data/annotation/Camarhynchus_parvulus_V1.1_genes.csv") 
  
get_genes <- function(comp.bed){
  parv.df.an.gr <- GRanges(parv.df.an$seqid, IRanges(as.numeric(parv.df.an$start), as.numeric(parv.df.an$end)))
  
  comp.bed <- as.data.frame(comp.bed)
  comp.bed.gr <- GRanges(comp.bed$CHROM, IRanges(comp.bed$start, comp.bed$end), peak = comp.bed$peak)
  
  comp_overlap <- findOverlaps(comp.bed.gr, parv.df.an.gr)
  
  comp.annotated <- cbind(comp.bed[queryHits(comp_overlap), ], parv.df.an[subjectHits(comp_overlap), c(-1, -2, -3)])
  
  comp.annotated <- comp.annotated %>% mutate(gene_name = gene_symbol) %>% 
    select(-gene_symbol, -Name)
  
  comp.annotated
}




