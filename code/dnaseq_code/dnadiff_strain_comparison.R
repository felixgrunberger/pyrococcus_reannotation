# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# compares two strains using output of Mummers dnadiff command
# Illumina/PacBio vs COM1 vs DSMZ
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(here)
library(data.table)
library(genoPlotR)
library(viridis)
library(ape)

# set root 
here()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load output from Mummer dnadiff
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
dsmznew_denovo  <- here("data/dnadiff_data/dsmznew_polished_dnadiff.mcoords")
com1_dsm3638    <- here("data/dnadiff_data/com1_dsm3638.dnadiff.mcoords")
dsm3638_dsmznew <- here("data/dnadiff_data/dsm3638_dsmznew.dnadiff.mcoords")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# preparation for plotting 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
com1_dsm3638_com1 <- fread(com1_dsm3638) %>%
  as.data.table() %>%
  dplyr:: mutate(name = V12,
                 start = V1,
                 end = V2,
                 strand = 1,
                 length = V2- V1,
                 pid = NA,
                 gene = NA,
                 synonym = NA,
                 product =NA,
                 proteinid = NA,
                 feature = NA,
                 gene_type = "side_blocks",
                 col = "white",
                 fill = viridis_pal(option = "D",begin = 0,end = 1)(10)[9], 
                 lty = 1,
                 lwd = .7,
                 pch = 20,
                 cex = -2) 
com1_dsm3638_com1 <- com1_dsm3638_com1[,-(1:13)]

com1_dsm3638_dsm3638 <- fread(com1_dsm3638) %>%
  as.data.table() %>%
  dplyr:: mutate(name = V12,
                 start = V3,
                 end = V4,
                 strand = 1,
                 length = V4- V3,
                 pid = NA,
                 gene = NA,
                 synonym = NA,
                 product =NA,
                 proteinid = NA,
                 feature = NA,
                 gene_type = "side_blocks",
                 col = "white",
                 fill = viridis_pal(option = "D",begin = 0,end = 1)(10)[6], 
                 lty = 1,
                 lwd = .5,
                 pch = 16,
                 cex = 3) 
com1_dsm3638_dsm3638 <- com1_dsm3638_dsm3638[,-(1:13)]

dsm3638_dsmznew_dsm3638 <- fread(dsm3638_dsmznew) %>%
  as.data.table() %>%
  dplyr:: mutate(name = V12,
                 start = V1,
                 end = V2,
                 strand = 1,
                 length = V2- V1,
                 pid = NA,
                 gene = NA,
                 synonym = NA,
                 product =NA,
                 proteinid = NA,
                 feature = NA,
                 gene_type = "side_blocks",
                 col = "white",
                 fill = viridis_pal(option = "D",begin = 0,end = 1)(10)[6], 
                 lty = 1,
                 lwd = .7,
                 pch = 20,
                 cex = -2) 
dsm3638_dsmznew_dsm3638 <- dsm3638_dsmznew_dsm3638[,-(1:13)]

dsm3638_dsmznew_dsmznew <- fread(dsm3638_dsmznew) %>%
  as.data.table() %>%
  dplyr:: mutate(name = V12,
                 start = V3,
                 end = V4,
                 strand = 1,
                 length = V4- V3,
                 pid = NA,
                 gene = NA,
                 synonym = NA,
                 product =NA,
                 proteinid = NA,
                 feature = NA,
                 gene_type = "side_blocks",
                 col = "white",
                 fill = viridis_pal(option = "D",begin = 0,end = 1)(10)[4], 
                 lty = 1,
                 lwd = .5,
                 pch = 16,
                 cex = 3) 
dsm3638_dsmznew_dsmznew <- dsm3638_dsmznew_dsmznew[,-(1:13)]

com1_vs_dsm3638 <- fread(com1_dsm3638) %>%
  as.data.table() %>%
  dplyr:: mutate(start1 = V1,
                 end1 = V2,
                 start2 = V3,
                 end2 = V4,
                 direction = 1,
                 col = "blue",
                 lwd = 10,
                 length = V2 - V1)
com1_vs_dsm3638 <- com1_vs_dsm3638[,-(1:13)]

dsm3638_vs_dsmznew <- fread(dsm3638_dsmznew) %>%
  as.data.table() %>%
  dplyr:: mutate(start1 = V1,
                 end1 = V2,
                 start2 = V3,
                 end2 = V4,
                 direction = 1,
                 col = "blue",
                 lwd = 10,
                 length = V2 - V1)
dsm3638_vs_dsmznew <- dsm3638_vs_dsmznew[,-(1:13)]

set_new <- c(list(as.comparison(com1_vs_dsm3638),as.comparison(dsm3638_vs_dsmznew)))
set_new[[1]]$col <- apply_color_scheme(set_new[[1]]$direction,
                                       color_scheme="grey",
                                       transparency=1)
set_new[[2]]$col <- apply_color_scheme(set_new[[2]]$direction,
                                       color_scheme="grey",
                                       transparency=1)

dna_set <- list(as.dna_seg(com1_dsm3638_com1),as.dna_seg(com1_dsm3638_dsm3638),as.dna_seg(dsm3638_dsmznew_dsmznew))
dna_set[[1]]$fill <- viridis_pal(option = "D",begin = 0,end = 1)(10)[9] 
dna_set[[2]]$fill <- viridis_pal(option = "D",begin = 0,end = 1)(10)[6] 
dna_set[[3]]$fill <- viridis_pal(option = "D",begin = 0,end = 1)(10)[4] 


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot with genoplotR (and save as PDF)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
pdf(here("figures/dnaseq_figures/strain_comparison.pdf"), 
    width = 6, height = 2, paper = "special",onefile=FALSE)
print(plot_gene_map(dna_segs=dna_set,
                    dna_seg_labels=c("com1","dsmz3638","new"),
                    comparisons = set_new,
                    dna_seg_scale=F, scale=T))
dev.off()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot position of IS elements
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
is_pfu_new <- read.gff(here("data/ises_data/pfu_dsmz_assembly_ises.fasta.gff")) %>%
  dplyr::filter(type == "insertion_sequence") %>%
  mutate(chr = "CP023154",
         value = 1,
         start_read = start,
         end_read = end,
         id = chr,
         width = end_read-start_read) %>%
  dplyr::select(chr, start_read, end_read, value, id, width)

is_pfu_old <- read.gff(here("data/ises_data/pfu_dsm3638_ises.fasta.gff")) %>%
  dplyr::filter(type == "insertion_sequence") %>%
  mutate(chr = "CP023154",
         value = 1,
         start_read = start,
         end_read = end,
         id = chr,
         width = end_read-start_read) %>%
  dplyr::select(chr, start_read, end_read, value, id, width)

pdf(here("figures/dnaseq_figures/is_position_dsm3638_old.pdf"), 
    width = 24, height = 3, paper = "special",onefile=FALSE)
ggplot(data = is_pfu_old, aes(y = 0)) +
  geom_segment(aes(x = start_read, xend = end_read, yend = 0), size = 30) +
  scale_x_continuous(limits=c(0,1908256)) +
  xlab("genomic position") +
  theme_minimal() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggtitle("dsm3638_old")
dev.off() 

pdf(here("figures/dnaseq_figures/is_position_dsm3638_new.pdf"), 
    width = 24, height = 3, paper = "special",onefile=FALSE)
ggplot(data = is_pfu_new, aes(y = 0)) +
  geom_segment(aes(x = start_read, xend = end_read, yend = 0), size = 30) +
  scale_x_continuous(limits=c(0,1889914)) +
  xlab("genomic position") +
  theme_minimal() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggtitle("dsm3638_new")
dev.off()

