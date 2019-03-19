# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Bidirectional transcripts, coverage density plots
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(here)
library(data.table)
library(viridis)
library(CoverageView)
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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load data from TSSpredator and genome gff file and make .bed files for further calculations
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> split dataset according to gene orientation 
split_primary <- function(type, selected_strand) {
  # filter for primary TSS
  tss_primary <- fread(here("data/annogesic_data/MasterTable.tsv")) %>%
    as.data.table() %>%
    filter(Primary == 1) %>%
    select(Pos, Strand, Locus_tag)
  
  # load gff annotation file and add information about divergent/non divergent orientation
  gff_table <- readGFF(file = here("data/genome_data/CP023154.gff")) %>% 
    as.data.frame() %>%
    dplyr::filter(type == "exon") %>%
    mutate(id = substr(Parent, 1,15)) %>%
    select(strand, id) %>%
    filter(isUnique(id)) %>%
    mutate(orientation = ifelse(strand == "+" & lag(strand) == "+", "non_divergent",
                                ifelse(strand == "-" & lead(strand) == "-", "non_divergent", 
                                       ifelse(strand == "+" & lag(strand) == "-", "divergent",
                                              ifelse(strand == "-" & lead(strand) == "+", "divergent", NA))))) 
  
  # join 2 tables and filter for orientation and strand
  left_join(tss_primary, gff_table, by = c("Locus_tag" = "id")) %>%
    dplyr::filter(orientation == type, strand == selected_strand) %>%
    mutate(chr = "CP023154",
           start = Pos,
           end = Pos,
           name = Locus_tag) %>%
    dplyr::select(chr,start,end,name)
}

# >> save files as .bed file
write.table(split_primary(type = "non_divergent", selected_strand = "+"),
            file = here("data/annogesic_data/non_divergent_plus.bed"), row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(split_primary(type = "non_divergent", selected_strand = "-"),
            file = here("data/annogesic_data/non_divergent_minus.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(split_primary(type = "divergent", selected_strand = "+"),
            file = here("data/annogesic_data/divergent_plus.bed"), row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(split_primary(type = "divergent", selected_strand = "-"),
            file = here("data/annogesic_data/divergent_minus.bed"),row.names = F, col.names = F, quote = F, sep = "\t")     


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# coverage plots using CoverageView
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> calculate coverage for sense/antisense strands for all sequencing data sets split by gene orientation
tss_coverage_calculation <- function(seq_type, divergent_type){
  # input sequencing bigwig files
  inputfile_f <- paste(here("data/annogesic_data/input/"),seq_type,"_in_CP023154_plus.bw", sep = "")
  inputfile_r <- paste(here("data/annogesic_data/input/"),seq_type,"_in_CP023154_minus.bw", sep = "")
  trm_f <- CoverageBigWigFile(inputfile_f)
  trm_r <- CoverageBigWigFile(inputfile_r)
  
  # calculate count matrix (dependet on divergent / non divergent type) input from previous function (1)
  # caluculate 400 bps up/downstream in a 10 bp window
  pfu_cov_f_sense    <- t(cov.matrix(trm_f,coordfile=paste(here("data/annogesic_data/"),divergent_type,"_plus.bed", sep = ""),
                                     extend = 400, num_cores = 4, bin_width = 10))
  pfu_cov_f_anti     <- t(cov.matrix(trm_r,coordfile=paste(here("data/annogesic_data/"),divergent_type,"_plus.bed", sep = ""),
                                     extend = 400, num_cores = 4, bin_width = 10))
  pfu_cov_r_anti     <- t(-cov.matrix(trm_f,coordfile=paste(here("data/annogesic_data/"),divergent_type,"_minus.bed", sep = ""),
                                      extend = 400, num_cores = 4, bin_width = 10))
  pfu_cov_r_sense    <- t(-cov.matrix(trm_r,coordfile=paste(here("data/annogesic_data/"),divergent_type,"_minus.bed", sep = ""),
                                      extend = 400, num_cores = 4, bin_width = 10))
  
  # bind values for both strands
  pfu_cov_sense <- rbind(pfu_cov_f_sense, pfu_cov_r_sense[,c(ncol(pfu_cov_r_sense):1)])
  pfu_cov_anti  <- rbind(pfu_cov_f_anti, pfu_cov_r_anti[,c(ncol(pfu_cov_r_anti):1)])
  
  # scale values by dividing each value/sum(values) for one TSS
  freq_sense <- scale(t(pfu_cov_sense), center = FALSE, 
                      scale = colSums(t(pfu_cov_sense)))
  freq_anti <- scale(t(pfu_cov_anti), center = FALSE, 
                     scale = colSums(t(pfu_cov_anti)))
  
  # make count matrix
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
frag_divergent    <- tss_coverage_calculation(seq_type = "FRAG", divergent_type = "divergent")
tex_divergent     <- tss_coverage_calculation(seq_type = "plus_TEX", divergent_type = "divergent")
notex_divergent   <- tss_coverage_calculation(seq_type = "minus_TEX", divergent_type = "divergent")

frag_non_divergent    <- tss_coverage_calculation(seq_type = "FRAG", divergent_type = "non_divergent")
tex_non_divergent     <- tss_coverage_calculation(seq_type = "plus_TEX", divergent_type = "non_divergent")
notex_non_divergent   <- tss_coverage_calculation(seq_type = "minus_TEX", divergent_type = "non_divergent")


divergent_coverage <- rbind(frag_divergent, tex_divergent, notex_divergent) %>%
  mutate(number = as.numeric(number),
         SENSE = as.numeric(SENSE),
         ANTISENSE = as.numeric(ANTISENSE))

non_divergent_coverage <- rbind(frag_non_divergent, tex_non_divergent, notex_non_divergent) %>%
  mutate(number = as.numeric(number),
         SENSE = as.numeric(SENSE),
         ANTISENSE = as.numeric(ANTISENSE))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot coverage plots (divergent & non-divergent)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> divergent
divergent_plot <- ggplot() + 
  stat_smooth(data = divergent_coverage, aes(x = number, y= SENSE, fill = TYPE), geom = "area", method = "loess", span = 0.07) +
  stat_smooth(data = divergent_coverage, aes(x = number, y= SENSE, color = TYPE, linetype= TYPE), geom = "line", method = "loess", span = 0.07, size = 2) +
  stat_smooth(data = divergent_coverage, aes(x = number, y= ANTISENSE, fill = TYPE), geom = "area", method = "loess", span = 0.07) +
  stat_smooth(data = divergent_coverage, aes(x = number, y= ANTISENSE, color = TYPE, linetype= TYPE), geom = "line", method = "loess", span = 0.07, size = 1.5) +
  theme_Publication() + 
  scale_linetype_manual(values=c("solid","twodash", "dotted")) +
  scale_fill_viridis(discrete = T, begin = 0.1, end = 0.7, alpha = 0.3) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.7, alpha = 1) +
  labs(x = "Position to pTSS [nt]", y = "average coverage") + 
  geom_vline(xintercept = c(0), linetype = "longdash", alpha = 0.5) +
  geom_vline(xintercept = -50, color = "red1", linetype = "longdash", alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05,0.12)) +
  ggtitle("head-to-head orientation") 

pdf(here("figures/rnaseq_figures/bidirectional_headtohead_cov.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
divergent_plot
dev.off()

# >> non_divergent
non_divergent_plot <- ggplot() + 
  stat_smooth(data = non_divergent_coverage, aes(x = number, y= SENSE, fill = TYPE), geom = "area", method = "loess", span = 0.07) +
  stat_smooth(data = non_divergent_coverage, aes(x = number, y= SENSE, color = TYPE, linetype= TYPE), geom = "line", method = "loess", span = 0.07, size = 2) +
  stat_smooth(data = non_divergent_coverage, aes(x = number, y= ANTISENSE, fill = TYPE), geom = "area", method = "loess", span = 0.07) +
  stat_smooth(data = non_divergent_coverage, aes(x = number, y= ANTISENSE, color = TYPE, linetype = TYPE), geom = "line", method = "loess", span = 0.07, size = 1.5) +
  theme_Publication() + 
  scale_linetype_manual(values=c("solid","twodash", "dotted")) +
  scale_fill_viridis(discrete = T, begin = 0.1, end = 0.7, alpha = 0.3) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.7, alpha = 1) +
  labs(x = "Position to pTSS [nt]", y = "average coverage") + 
  geom_vline(xintercept = c(0), linetype = "longdash", alpha = 0.5) +
  geom_vline(xintercept = -50, color = "red1", linetype = "longdash", alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05,0.12)) +
  ggtitle("head-to-tail orientation") 

pdf(here("figures/rnaseq_figures/bidirectional_headtotail_cov_guide.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
non_divergent_plot
dev.off()


