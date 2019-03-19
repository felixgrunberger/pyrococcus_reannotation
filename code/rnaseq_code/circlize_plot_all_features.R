# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Circlize plot showing all features from DNA and RNA seq
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(viridis)
library(circlize)
library(Biostrings)
library(data.table)
library(ape)
library(stringr)
library(readr)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Prepare data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> GENOME INFORMATION
dsm3638_gff = here("data/genome_data/CP023154.gff")
dsm3638_gff_file <- read.gff(dsm3638_gff)
gff_table <- read.table(file = dsm3638_gff, 
                        sep = "\t", header = FALSE, stringsAsFactors = F, fill = F, quote = "", 
                        col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
  dplyr::filter(type == "gene") %>%
  mutate(chr = "CP023154",
         start = start,
         end = end,
         id = substr(attributes, 4,18)) %>%
  dplyr::select(chr, start, end, id)

gff_table_plus <- read.table(file = dsm3638_gff, 
                             sep = "\t", header = FALSE, stringsAsFactors = F, fill = F, quote = "", 
                             col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
  dplyr::filter(type == "gene" & strand == "+") %>%
  mutate(chr = "CP023154",
         start = start,
         end = end,
         id = substr(attributes, 4,18)) %>%
  dplyr::select(chr, start, end, id)

gff_table_minus <- read.table(file = dsm3638_gff, 
                              sep = "\t", header = FALSE, stringsAsFactors = F, fill = F, quote = "", 
                              col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
  dplyr::filter(type == "gene" & strand == "-") %>%
  mutate(chr = "CP023154",
         start = start,
         end = end,
         id = substr(attributes, 4,18)) %>%
  dplyr::select(chr, start, end, id)

gff_table_strands <- list(gff_table_plus, gff_table_minus)

gff_table_distance <- read.table(file = dsm3638_gff, 
                                 sep = "\t", header = FALSE, stringsAsFactors = F, fill = F, quote = "", 
                                 col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
  dplyr::filter(type == "gene") %>%
  mutate(chr = "CP023154",
         start = start,
         end = end,
         start_next = c(start[-1],NA),
         dist = start_next-end) %>%
  dplyr::select(chr, start, end, dist)


# mastertable tsspredator
mastertable <- fread(here("data/annogesic_data/CP023154_merge_features.gff")) %>%
  filter(V3 != "contig")

# ----------------------- #
########################### DNA based annotation
# 1. gc content (histogram)
# ------------------------------------------------------------------------------------------------------------------------- #
# 1. GC VALUE
dsm3639_fa <- here("data/fasta_data/pfu_dsmz_assembly.fasta")
fasta <- readDNAStringSet(file = dsm3639_fa)
names(fasta) <- "CP023154"
stepwidth <- 1000
gc_value <- matrix(ncol = 4, nrow = sum(1:length(fasta$CP023154) %% stepwidth == 0)) %>%
  as.data.table() %>%
  mutate(chr   = "CP023154",
         start =  which(1:length(fasta$CP023154) %% stepwidth == 0) - stepwidth + 1,
         end   =  which(1:length(fasta$CP023154) %% stepwidth == 0) + 1,
         value    = 0) %>%
  dplyr::select(chr, start, end, value)

for ( i in 1:(dim(gc_value)[1])){
  gc_value$value[i] <-  GC.content(as.DNAbin(fasta$CP023154[gc_value$start[i]:gc_value$end[i]]))
}

gc_max_height = max(gc_value$value)  
gc_min_height =  min(gc_value$value)
gc_mean_height <- mean(gc_value$value)
at_value <- gc_value %>%
  mutate(value = 1-value)

at_max_height = max(at_value$value)  
at_min_height =  min(at_value$value)

gc_value$value1 <- at_value$value
gc_value$value1 <- NULL

base_values <- list(gc_value, at_value)
gc_value$value_1 <- NA
gc_value$value_2 <- NA

for (i in 1:length(gc_value$value)){
  if (gc_value$value[i] >= mean(gc_value$value)){
    gc_value$value_1[i] <- gc_value$value[i]
    gc_value$value_2[i] <- NA
  } 
  if (gc_value$value[i] < mean(gc_value$value)){
    gc_value$value_1[i] <- NA
    gc_value$value_2[i] <- gc_value$value[i]
  } 
}
gc_value2 <- gc_value
gc_value2$value <- NULL
gc_value3 <- gc_value2
gc_value2$value_2 <- NULL
gc_value3$value_1 <- NULL
gc_values_all <- list(gc_value2, gc_value3)

# 2. CDS plus/minus (rect-plot)
# 3. genomic density (rainfall plot)
# 4. CRISPR
CP023154_crispr <-  dsm3638_gff_file %>%
  dplyr::filter(type == "repeat_region") %>%
  mutate(chr = "CP023154",
         start = start,
         end = end,
         id = substr(attributes, 4,18)) %>%
  dplyr::select(chr, start, end, id)

# 5. tRNA
CP023154_tRNA <-  dsm3638_gff_file %>%
  dplyr::filter(type == "tRNA") %>%
  mutate(chr = "CP023154",
         start = start,
         end = end,
         id = substr(attributes, 4,18)) %>%
  dplyr::select(chr, start, end, id)

# 6. rRNA
CP023154_rRNA <-  dsm3638_gff_file %>%
  dplyr::filter(type == "rRNA") %>%
  mutate(chr = "CP023154",
         start = start,
         end = end,
         id = substr(attributes, 4,18)) %>%
  dplyr::select(chr, start, end, id)

CP023154_levels <- list(CP023154_crispr,
                        CP023154_tRNA, 
                        CP023154_rRNA)

# ----------------------- #
########################### RNA based annotation
# 6. TSS
CP023154_TSS <- mastertable %>%
  dplyr::filter(V3 == "TSS") %>%
  mutate(chr = "CP023154",
         start = V4,
         end = V5,
         id = V2) %>%
  dplyr::select(chr, start, end, id)

# 7. 5UTR 
CP023154_5UTR <- mastertable %>%
  dplyr::filter(V3 == "5UTR") %>%
  mutate(chr = "CP023154",
         start = V4,
         end = V5,
         value = as.numeric(str_split_fixed(str_split_fixed(V9,";",8)[,4], "=",2)[,2])
  ) %>%
  dplyr::select(chr, start, end, value)

# 8. 3UTR  
CP023154_3UTR <- mastertable %>%
  dplyr::filter(V3 == "3UTR") %>%
  mutate(chr = "CP023154",
         start = V4,
         end = V5,
         value = as.numeric(str_split_fixed(str_split_fixed(V9,";",6)[,3], "=",2)[,2])
  ) %>%
  dplyr::select(chr, start, end, value)

# 7. circRNA
CP023154_circRNA <- mastertable %>%
  dplyr::filter(V3 == "circRNA") %>%
  mutate(chr = "CP023154",
         start = V4,
         end = V5,
         id = V2) %>%
  dplyr::select(chr, start, end, id)

# 8. ncRNA
CP023154_ncRNA <- mastertable %>%
  dplyr::filter(V3 == "ncRNA") %>%
  mutate(chr = "CP023154",
         start = V4,
         end = V5,
         id = V2) %>%
  dplyr::select(chr, start, end, id)

# 9. operon
CP023154_operon <- mastertable %>%
  dplyr::filter(V3 == "operon") %>%
  mutate(chr = "CP023154",
         start = V4,
         end = V5,
         id = V2) %>%
  dplyr::select(chr, start, end, id)

# 10. sORF
CP023154_sORF <- mastertable %>%
  dplyr::filter(V3 == "sORF") %>%
  mutate(chr = "CP023154",
         start = V4,
         end = V5,
         id = V2) %>%
  dplyr::select(chr, start, end, id)

# 11. terminator
CP023154_terminator <- mastertable %>%
  dplyr::filter(V3 == "terminator") %>%
  mutate(chr = "CP023154",
         start = V4,
         end = V5,
         id = V2) %>%
  dplyr::select(chr, start, end, id)

CP023154_rna_levels <- list(CP023154_TSS, 
                            CP023154_circRNA, 
                            CP023154_ncRNA, 
                            CP023154_operon, 
                            CP023154_sORF, 
                            CP023154_terminator)

# 12. ISelements calculated by ISEScan-1.6
ises <- read.gff(here("data/ises_data/pfu_dsmz_assembly_ises.fasta.gff")) %>%
  dplyr::filter(type == "insertion_sequence") %>%
  mutate(chr = "CP023154",
         value = 1) %>%
  dplyr::select(chr, start, end, value)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
dev.off()

pdf(here("figures/rnaseq_figures/circlize_90.pdf"), 
    width = 10, height = 10, paper = "special",onefile=FALSE)

# ---- #
# DELETE PREVIOUS PLOT
circos.clear()

# ---- #
# SET LINE WIDTH RELAITVE TO DEFAULT
par(lwd = 0.5)

# ---- #
# ARRANGE CIRCULAR LAYOUT: with padding of 0 to all sites and a gap degree of 45 degrees, start degree
circos.par("cell.padding" = c(0, 0, 0, 0), "gap.degree" = 90, "start.degree" = 0)

# ---- #
# Initialize the circular layout with genomic data
circos.genomicInitialize(gff_table)

# ---------------------------------------------------------------------------------------------------------------- #
# 1. GC VALUE
color_gc_gradient <- rev(viridis_pal(option = "D")(24)[c(12,2)])
circos.genomicTrackPlotRegion(gc_values_all, ylim = c(gc_min_height, gc_max_height), panel.fun = function(region, value, ...) {
  i = getI(...)
  circos.genomicLines(region, value, 
                      area = F, lwd = .9, baseline = mean(gc_value$value),
                      col = color_gc_gradient[i], type = "h", border = NA, ...)
}, bg.border = NA, track.height = 0.10)

# -------- #
# 2 random colorsets of the virdis palette --> 1 lighter and 1 darker
rand_col = function(k) {
  return(viridis_pal(option = "D")(24)[4:8][1 +3*runif(k)])
}
rand_col2 = function(k) {
  return(viridis_pal(option = "D")(24)[8:14][1+3*runif(k)])
}

# ---------------------------------------------------------------------------------------------------------------- #
# 2. genomicrect plot of genes annotated on both strand
circos.genomicTrackPlotRegion(gff_table_strands, stack = TRUE, panel.fun = function(region, value, ...){
  i <- getI(...)
  if (i == 1) {
    circos.genomicRect(region, value, col = rand_col2(nrow(region)), border = NA, ...)
  } else {
    circos.genomicRect(region, value, col = rand_col(nrow(region)), border = NA, ...)
  }
}, bg.border = NA, track.height = 0.05, track.margin = c(0,0.02))

# ---------------------------------------------------------------------------------------------------------------- #
# 3. GENOMIC RAINFALL (--> density)
circos.genomicRainfall(gff_table_distance, stack = F, ylim = c(-2.5,2.5), track.height = 0.10, track.margin =c(0,0.05),
                       pch = 16, bg.border = F, cex = 1, col = viridis_pal(alpha = 0.20,option="D")(24)[10])

# ---------------------------------------------------------------------------------------------------------------- #
# 4. CRISPR | 5. tRNA | 6. rRNA
# different levels of dsm3638 annotation (direct_repeat, pseudogene, rRNA, snoRNA, SRP_RNA, tRNA, TSS, ncRNA, terminator)
color_rect2 = rev(viridis_pal(option = "D")(24)[c(2,8,12)])
circos.genomicTrackPlotRegion(CP023154_levels, stack = TRUE,panel.fun = function(region, value, ...) {
  i = getI(...)
  circos.genomicRect(region, value, col = NA, lwd = 1,border = color_rect2[i], ...)
}, bg.border = "grey",bg.lwd=1,track.height = 0.10, track.margin = c(0,0))


# ---------------------------------------------------------------------------------------------------------------- #
# 7. circRNA 8. ncRNA 9. operon 10. sORF 11. terminator
color_rect_tss = rev(viridis_pal(option = "D")(24)[c(3,5,7,9,11,13,15)])
circos.genomicTrackPlotRegion(CP023154_rna_levels, stack = TRUE,panel.fun = function(region, value, ...) {
  i = getI(...)
  circos.genomicRect(region, value, col = NA, lwd = 0.5,border = color_rect_tss[i], ...)
},track.height = 0.2, track.margin = c(0,0))

# ---------------------------------------------------------------------------------------------------------------- #
# 12. ises
circos.genomicTrackPlotRegion(ises, stack = TRUE,
                              panel.fun = function(region, value, ...){
                                circos.genomicRect(region, value, col = "black", lwd = 0.5, ...)
                              },track.height = 0.1, track.margin = c(0,0))

# ----------------- #
dev.off()
circos.clear()