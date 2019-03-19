# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Average read coverage of 3 data sets in relative position to 4 TSS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(CoverageView)
library(data.table)
library(tidyverse)
library(here)
library(stringr)
library(viridis)
library(Biobase)
library(ggthemes)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# theme for plotting
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
theme_Publication <- function(base_size=14) {
  (theme_foundation(base_size=base_size, base_family="Helvetica")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> strand specific wig files (input for ANNOgesic) were converted using wigToBigWig (https://anaconda.org/bioconda/ucsc-wigtobigwig)

# >> load TSS and write to bed file for each category
get_tss_bed <- function(what_type, strand){
  return(fread(input = here("data/annogesic_data/CP023154_TSS.gff"),skip = 1) %>%
           as.data.table() %>%
           dplyr::filter(V7 == strand) %>%
           mutate(type = str_split(V9, ";",n = 7, simplify = T)[,3]) %>%
           dplyr::filter(type == paste("type=",what_type, sep="")) %>%
           dplyr::mutate(chr   = V1,
                         start = V4,
                         end   = V5,
                         name = type) %>%
           dplyr::select(chr,start,end,name))
}

# >> write positions to table
write.table(get_tss_bed("Primary", "+"),  file = here("data/annogesic_data/CP023154_TSS_Primary_p.gff.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(get_tss_bed("Secondary","+"), file = here("data/annogesic_data/CP023154_TSS_Secondary_p.gff.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(get_tss_bed("Internal", "+"), file = here("data/annogesic_data/CP023154_TSS_Internal_p.gff.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(get_tss_bed("Antisense", "+"),file = here("data/annogesic_data/CP023154_TSS_Antisense_p.gff.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(get_tss_bed("Primary", "-"),  file = here("data/annogesic_data/CP023154_TSS_Primary_m.gff.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(get_tss_bed("Secondary","-"), file = here("data/annogesic_data/CP023154_TSS_Secondary_m.gff.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(get_tss_bed("Internal", "-"), file = here("data/annogesic_data/CP023154_TSS_Internal_m.gff.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(get_tss_bed("Antisense", "-"),file = here("data/annogesic_data/CP023154_TSS_Antisense_m.gff.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# calculate average read density for each TSS category
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
tss_coverage <- function(seq_type, tss_type, filtered = F){
  # load bw file for each sequencing type
  inputfile_f <- paste(here("data/annogesic_data/input/"),seq_type,"_in_CP023154_plus.bw", sep = "")
  inputfile_r <- paste(here("data/annogesic_data/input/"),seq_type,"_in_CP023154_minus.bw", sep = "")
  # use coverageBigWig to calculate coverage
  trm_f <- CoverageBigWigFile(inputfile_f)
  trm_r <- CoverageBigWigFile(inputfile_r)
  # make a strand specific count matrix
  pfu_cov_f_sense    <- t(cov.matrix(trm_f,coordfile=paste(here("data/annogesic_data/CP023154_TSS_"),tss_type,"_p.gff.bed",sep = ""),extend=400,num_cores=2, bin_width=10))
  pfu_cov_f_anti     <- t(cov.matrix(trm_r,coordfile=paste(here("data/annogesic_data/CP023154_TSS_"),tss_type,"_p.gff.bed",sep = ""),extend=400,num_cores=2, bin_width=10))
  pfu_cov_r_anti     <- t(-cov.matrix(trm_f,coordfile=paste(here("data/annogesic_data/CP023154_TSS_"),tss_type,"_m.gff.bed",sep = ""),extend=400,num_cores=2, bin_width=10))
  pfu_cov_r_sense    <- t(-cov.matrix(trm_r,coordfile=paste(here("data/annogesic_data/CP023154_TSS_"),tss_type,"_m.gff.bed",sep = ""),extend=400,num_cores=2, bin_width=10))
  # combine values for both strands
  pfu_cov_sense <- rbind(pfu_cov_f_sense, pfu_cov_r_sense[,c(ncol(pfu_cov_r_sense):1)])
  pfu_cov_anti  <- rbind(pfu_cov_f_anti, pfu_cov_r_anti[,c(ncol(pfu_cov_r_anti):1)])
  # scale values by dividing each value/sum(values) 
  freq_sense <- scale(t(pfu_cov_sense), center = FALSE, 
                      scale = colSums(t(pfu_cov_sense)))
  freq_anti <- scale(t(pfu_cov_anti), center = FALSE, 
                     scale = colSums(t(pfu_cov_anti)))
  # counts were calculated in a window from -400 to +400 to a TSS in a 10 bp window
  line <- seq(from = -400+5, to = 400-5, by = 10)
  sense <- rowMeans(freq_sense)
  anti <- -rowMeans(freq_anti)
  data_TLE <- matrix(ncol = 4, nrow = 400/10*2)
  data_TLE[,1] <- line
  data_TLE[,2] <- sense
  data_TLE[,3] <- anti
  data_TLE[,4] <- seq_type
  data_TLE <- as.data.table(data_TLE)
  colnames(data_TLE) <- c("number", "SENSE", "ANTISENSE", "TYPE")
  return(data_TLE)
}

# >> calculate coverage & format for ggploting, primary
frag_cov_p <- tss_coverage(seq_type = "FRAG", tss_type = "Primary")
tex_cov_p  <- tss_coverage(seq_type = "plus_TEX", tss_type = "Primary")
notex_cov_p  <- tss_coverage(seq_type = "minus_TEX", tss_type = "Primary")
# >> calculate coverage & format for ggploting, secondary
frag_cov_s <- tss_coverage(seq_type = "FRAG", tss_type = "Secondary")
tex_cov_s  <- tss_coverage(seq_type = "plus_TEX", tss_type = "Secondary")
notex_cov_s  <- tss_coverage(seq_type = "minus_TEX", tss_type = "Secondary")
# >> calculate coverage & format for ggploting, internal
frag_cov_i <- tss_coverage(seq_type = "FRAG", tss_type = "Internal")
tex_cov_i  <- tss_coverage(seq_type = "plus_TEX", tss_type = "Internal")
notex_cov_i  <- tss_coverage(seq_type = "minus_TEX", tss_type = "Internal")
# >> calculate coverage & format for ggploting, antisense 
frag_cov_a <- tss_coverage(seq_type = "FRAG", tss_type = "Antisense")
tex_cov_a  <- tss_coverage(seq_type = "plus_TEX", tss_type = "Antisense")
notex_cov_a  <- tss_coverage(seq_type = "minus_TEX", tss_type = "Antisense")

# >> combine three sequencing data sets for each frame
cov_p <- rbind(frag_cov_p, tex_cov_p, notex_cov_p) %>%
  mutate(number = as.numeric(number),
         SENSE = as.numeric(SENSE),
         ANTISENSE = as.numeric(ANTISENSE))
cov_s <- rbind(frag_cov_s, tex_cov_s, notex_cov_s) %>%
  mutate(number = as.numeric(number),
         SENSE = as.numeric(SENSE),
         ANTISENSE = as.numeric(ANTISENSE))
cov_i <- rbind(frag_cov_i, tex_cov_i, notex_cov_i) %>%
  mutate(number = as.numeric(number),
         SENSE = as.numeric(SENSE),
         ANTISENSE = as.numeric(ANTISENSE))
cov_a <- rbind(frag_cov_a, tex_cov_a, notex_cov_a) %>%
  mutate(number = as.numeric(number),
         SENSE = as.numeric(SENSE),
         ANTISENSE = as.numeric(ANTISENSE))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot coverage files
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> primary
pdf(here("figures/rnaseq_figures/coverage_3primary.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot() + 
  stat_smooth(data = cov_p, aes(x = number, y= SENSE, fill = TYPE), geom = "area", method = "loess", span = 0.07) +
  stat_smooth(data = cov_p, aes(x = number, y= SENSE, color = TYPE, linetype= TYPE), geom = "line", method = "loess", span = 0.07, size = 2) +
  theme_Publication() + 
  scale_linetype_manual(values=c("solid","twodash", "dotted")) +
  scale_fill_viridis(discrete = T, begin = 0.2, end = 0.6, alpha = 0.3) +
  scale_color_viridis(discrete = T, begin = 0.2, end = 0.6, alpha = 1) +
  labs(x = "Position to pTSS [nt]", y = "average coverage") + 
  geom_vline(xintercept = c(0), linetype = "longdash", alpha = 0.5) 
dev.off()

# >> secondary
pdf(here("figures/rnaseq_figures/coverage_3sec.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot() + 
  stat_smooth(data = cov_s, aes(x = number, y= SENSE, fill = TYPE), geom = "area", method = "loess", span = 0.07) +
  stat_smooth(data = cov_s, aes(x = number, y= SENSE, color = TYPE, linetype= TYPE), geom = "line", method = "loess", span = 0.07, size = 2) +
  theme_Publication() +  
  scale_linetype_manual(values=c("solid","twodash", "dotted")) +
  scale_fill_viridis(discrete = T, begin = 0.2, end = 0.6, alpha = 0.3) +
  scale_color_viridis(discrete = T, begin = 0.2, end = 0.6, alpha = 1) +
  labs(x = "Position to sTSS [nt]", y = "average coverage") + 
  geom_vline(xintercept = c(0), linetype = "longdash", alpha = 0.5) 
dev.off()

# >> internal
pdf(here("figures/rnaseq_figures/coverage_3internal.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot() + 
  stat_smooth(data = cov_i, aes(x = number, y= SENSE, fill = TYPE), geom = "area", method = "loess", span = 0.07) +
  stat_smooth(data = cov_i, aes(x = number, y= SENSE, color = TYPE, linetype= TYPE), geom = "line", method = "loess", span = 0.07, size = 2) +
  theme_Publication() + 
  scale_linetype_manual(values=c("solid","twodash", "dotted")) +
  scale_fill_viridis(discrete = T, begin = 0.2, end = 0.6, alpha = 0.3) +
  scale_color_viridis(discrete = T, begin = 0.2, end = 0.6, alpha = 1) +
  labs(x = "Position to iTSS [nt]", y = "average coverage") + 
  geom_vline(xintercept = c(0), linetype = "longdash", alpha = 0.5) 
dev.off()

# >> antisense
pdf(here("figures/rnaseq_figures/coverage_3antisense.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot() + 
  stat_smooth(data = cov_a, aes(x = number, y= SENSE, fill = TYPE), geom = "area", method = "loess", span = 0.07) +
  stat_smooth(data = cov_a, aes(x = number, y= SENSE, color = TYPE, linetype= TYPE), geom = "line", method = "loess", span = 0.07, size = 2) +
  theme_Publication() + 
  scale_linetype_manual(values=c("solid","twodash", "dotted")) +
  scale_fill_viridis(discrete = T, begin = 0.2, end = 0.6, alpha = 0.3) +
  scale_color_viridis(discrete = T, begin = 0.2, end = 0.6, alpha = 1) +
  labs(x = "Position to aTSS [nt]", y = "average coverage") + 
  geom_vline(xintercept = c(0), linetype = "longdash", alpha = 0.5) 
dev.off()
