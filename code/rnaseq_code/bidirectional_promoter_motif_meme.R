# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Bidirectional transcripts, promoter motifs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(here)
library(data.table)
library(viridis)
library(seqinr)
library(ggseqlogo)
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
# preapare data sets for analysis with meme
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
tss_coverage_orientation_anti <- function(seq_type, divergent_type){
  
  # input sequencing bigwig files
  inputfile_f <- paste(here("data/annogesic_data/input/"),seq_type,"_in_CP023154_plus.bw", sep = "")
  inputfile_r <- paste(here("data/annogesic_data/input/"),seq_type,"_in_CP023154_minus.bw", sep = "")
  trm_f <- CoverageBigWigFile(inputfile_f)
  trm_r <- CoverageBigWigFile(inputfile_r)
  
  # calculate coverage only of antisense files
  pfu_cov_f_anti     <- t(cov.matrix(trm_r,coordfile=paste("data/annogesic_data/",divergent_type,"_plus.bed", sep = ""),
                                     extend = 400, num_cores = 4, bin_width = 10))
  pfu_cov_r_anti     <- t(-cov.matrix(trm_f,coordfile=paste("data/annogesic_data/",divergent_type,"_minus.bed", sep = ""),
                                      extend = 400, num_cores = 4, bin_width = 10))
  
  # combine files
  pfu_cov_anti  <- rbind(pfu_cov_f_anti, pfu_cov_r_anti[,c(ncol(pfu_cov_r_anti):1)])
}

tex_non_divergent_anti <- tss_coverage_orientation_anti(seq_type = "plus_TEX", 
                                                        divergent_type = "non_divergent")
# scaled counts
tex_anti <- t(scale(t(tex_non_divergent_anti), center = FALSE, 
                    scale = colSums(t(tex_non_divergent_anti))))

# calculate sum and mean per gene (same values because of scaling in previous step)
tex_anti_frame <- tex_anti %>%
  as.tibble() %>%
  mutate(gene_sum = rowSums(.),
         gene_mean = rowMeans(.))

# filter for relevant position 100 bps upstream of pTSS, e.g. here 30:40
# calculate % of reads that can are found in this region (subset_sum)
# subtract expected % /11*gene_mean)
tex_anti_frame_filtered <- tex_anti_frame[,c(30:40,81,82)] %>%
  mutate(subset_sum = rowSums(.) - gene_sum + gene_mean,
         subset_diff = subset_sum - (11*gene_mean)) 

# combine with upstream sequence of pTSS
tss_gff_distance_non_divergent_plus <- tss_gff_distance %>%
  filter(!is.na(Pos)) %>%
  filter(orientation == "non_divergent", strand == "+")
tss_gff_distance_non_divergent_minus <- tss_gff_distance %>%
  filter(!is.na(Pos)) %>%
  filter(orientation == "non_divergent", strand == "-")

# read in fasta
pfu_fasta <- readDNAStringSet(file = here("data/genome_data/pfu_dsmz_assembly.fasta"))
names(pfu_fasta) <- "CP023154"
#
tss_gff_distance_non_divergent_all <- rbind(tss_gff_distance_non_divergent_plus,
                                            tss_gff_distance_non_divergent_minus)

tx_nd_anti <- cbind(tex_anti_frame_filtered,tss_gff_distance_non_divergent_all) %>%
  rowwise() %>%
  mutate(sequence = ifelse(strand == "+",as.character(pfu_fasta$CP023154[(Pos-50):(Pos)]),
                           ifelse(strand == "-", as.character(reverseComplement(pfu_fasta$CP023154[(Pos):(Pos+50)])))))

atss <- fread(here("data/annogesic_data/MasterTable.tsv")) %>%
  as.data.table() %>%
  filter(Antisense == 1) %>%
  mutate(antisense_start = Pos) %>%
  select(antisense_start, Strand, Locus_tag, Antisense)

tx_nd_anti_atss <- left_join(tx_nd_anti, atss, by = c("name" = "Locus_tag")) %>%
  rename("primary_start" = Pos)

# >> select names of diff > median
tx_nd_anti_atss_high <- tx_nd_anti_atss %>%
  dplyr::filter(subset_diff > 0.40)

write.fasta(as.list(tx_nd_anti_atss_high$sequence), 
            as.list(tx_nd_anti_atss_high$name), 
            file = here("data/meme_data/bidirectional.fasta")) 

# >> meme command (in terminal)
# meme "data/meme_data/bidirectional.fasta" \
# -oc "data/meme_data/bidrectional" \
# -dna -nostatus -time 18000 -mod zoops \
# -nmotifs 3 -minw 6 -maxw 50 -revcomp 



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load fasta from HTML output of `meme`
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
f_p <- read.table(here("data/meme_data/bidirectional_motif1.fasta"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1))
seq_p <- c(as.character(droplevels(f_p$V1)))
seq_p_bottom <- c(as.character(droplevels(reverseComplement(DNAStringSet(f_p$V1)))))
seq <- list(seq_p, seq_p_bottom)
names(seq) <- c("top", "bottom")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot motifs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> colors for bases
color_scale = make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                              cols=viridis_pal(begin = 0.2, end = 0.6,option = "viridis")(25)[c(1,7,13,19)])

# >> top strand
pdf(here("figures/rnaseq_figures/bidirectional_motif_top.pdf"), 
    width = 10, height = 3.5, paper = "special",onefile=FALSE)
ggplot() + 
  geom_logo(seq_p, font = "helvetica_regular", col_scheme = color_scale, seq_type = "dna") + 
  theme_logo() 
dev.off()

# >> bottom strand
pdf(here("figures/rnaseq_figures/bidirectional_motif_bottom.pdf"), 
    width = 10, height = 3.5, paper = "special",onefile=FALSE)
ggplot() + 
  geom_logo(seq_p_bottom, font = "helvetica_regular", col_scheme = color_scale, seq_type = "dna") + 
  theme_logo() 
dev.off()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# calculate position of motifs starting with output from meme
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
primary_pos_p <- fread(here("data/meme_data/bidirectional_positions.tsv"), fill = T) %>%
  filter(V2 == "+") %>%
  mutate(position = as.numeric(V3) - 50+11) %>% #(-50 of searched sequence + plus of motif width)
  select(position) 
primary_pos_n <- fread(here("data/meme_data/bidirectional_positions.tsv"), fill = T) %>%
  filter(V2 == "-") %>%
  mutate(position = as.numeric(V3) - 50+11) %>% #(-50 of searched sequence + plus normalized for position of third T)
  select(position) 
primary_pos <- rbind(primary_pos_p, primary_pos_n)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot motif positions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
pdf(here("figures/rnaseq_figures/bidirectional_histogram.pdf"), 
    width = 10, height = 2.5, paper = "special",onefile=FALSE)
ggplot(data = primary_pos, aes(x = position, y = ..count..)) +
  theme_Publication() +
  geom_histogram(aes(y=..count.., fill = ..count..),      # Histogram with density instead of count on y-axis
                 binwidth=1, size = 0.3, alpha = 0.6,colour = "white") +
  scale_fill_viridis(alpha = 0.3,discrete = F, option = "viridis", begin = 0.2, end = 0.6) +
  guides(fill = F) +
  ylab("number") +
  xlab("position to pTSS [nt]") +
  theme(text = element_text(size = 26))
dev.off()



