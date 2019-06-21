# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Enrichment of aTSS around ISelements
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(viridis)
library(Biostrings)
library(data.table)
library(ape)
library(ggthemes)
library(here)

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
# load data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> annotated aTSS
mastertable <- here("data/annogesic_data/MasterTable.tsv")
tss_table_antisense <- fread(mastertable) %>%
  as.data.table() %>%
  dplyr::filter(Antisense == 1) %>%
  mutate(type = "Antisense") %>%
  dplyr::select(Locus_tag, Strand, Pos,type, UTRlength) 


# >> basic genome file
gff_table <- read.table(file = here("data/genome_data/CP023154.gff"), 
                        sep = "\t", header = FALSE, stringsAsFactors = F, fill = F, quote = "",
                        col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
  dplyr::filter(type == "exon") %>%
  mutate(length = end - start,
         start = ifelse(start < end, start, end),
         chr = "CP023154",
         lag_strand = lag(strand),
         lead_strand = lead(strand),
         lag_end = lag(end),
         lead_start = lead(start),
         same_strand = ifelse(strand == "+" & strand == lag_strand, T, 
                              ifelse(strand == "-" & strand == lead_strand, T, F)),
         same_distance = ifelse(strand == "+" & strand == lag_strand, start - lag_end,
                                ifelse(strand =="-" & strand == lead_strand, lead_start - end,NA)),
         id = substr(attributes, 8,22)) %>%
  dplyr::select(id, start, end, strand, length, lag_strand, lead_strand, same_strand, same_distance)

# >> is element data from ISEScan1-6
is_pfu_new <- read.gff(here("data/ises_data/pfu_dsmz_assembly_ises.fasta.gff")) %>%
  dplyr::filter(type == "insertion_sequence") %>%
  mutate(chr = "CP023154",
         value = 1,
         start_read = start,
         end_read = end,
         id = chr) %>%
  dplyr::select(chr, start_read, end_read, value, id)

# >> find gene next to is element
is_pfu_new$cor_gene <- NA
for(i in 1:length(is_pfu_new$id)){
  if(which.min(abs(is_pfu_new$start_read[i]-gff_table$start)) == which.min(abs(is_pfu_new$start_read[i]-gff_table$end)) &&
     is_pfu_new$start_read[i] >= gff_table$start[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))] && is_pfu_new$start_read[i] <= gff_table$end[which.min(abs(is_pfu_new$start_read[i]-gff_table$end))]){
    is_pfu_new$cor_gene[i] <- gff_table$id[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))]
  }
  if(which.min(abs(is_pfu_new$start_read[i]-gff_table$start)) == which.min(abs(is_pfu_new$start_read[i]-gff_table$end)) &&
     is_pfu_new$start_read[i] >= gff_table$start[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))] && is_pfu_new$start_read[i] >= gff_table$end[which.min(abs(is_pfu_new$start_read[i]-gff_table$end))]){
    is_pfu_new$cor_gene[i] <- gff_table$id[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))]
  }
  if(which.min(abs(is_pfu_new$start_read[i]-gff_table$start)) == which.min(abs(is_pfu_new$start_read[i]-gff_table$end)) &&
     is_pfu_new$start_read[i] <= gff_table$start[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))] && is_pfu_new$start_read[i] <= gff_table$end[which.min(abs(is_pfu_new$start_read[i]-gff_table$end))]){
    is_pfu_new$cor_gene[i] <- gff_table$id[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))]
  }
  if(which.min(abs(is_pfu_new$start_read[i]-gff_table$start)) > which.min(abs(is_pfu_new$start_read[i]-gff_table$end)) &&
     is_pfu_new$start_read[i] <= gff_table$start[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))] && is_pfu_new$start_read[i] >= gff_table$end[which.min(abs(is_pfu_new$start_read[i]-gff_table$end))]){
    is_pfu_new$cor_gene[i] <- gff_table$id[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))]
  }
  if(which.min(abs(is_pfu_new$start_read[i]-gff_table$start)) > which.min(abs(is_pfu_new$start_read[i]-gff_table$end)) &&
     is_pfu_new$start_read[i] >= gff_table$start[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))] && is_pfu_new$start_read[i] <= gff_table$end[which.min(abs(is_pfu_new$start_read[i]-gff_table$end))+1]){
    is_pfu_new$cor_gene[i] <- gff_table$id[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))]
  }
  if(which.min(abs(is_pfu_new$start_read[i]-gff_table$start)) > which.min(abs(is_pfu_new$start_read[i]-gff_table$end)) &&
     is_pfu_new$start_read[i] < gff_table$start[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))] && is_pfu_new$start_read[i] < gff_table$end[which.min(abs(is_pfu_new$start_read[i]-gff_table$end))]){
    is_pfu_new$cor_gene[i] <- gff_table$id[which.min(abs(is_pfu_new$start_read[i]-gff_table$start))]
  }
}


# >> combine TSS and IS file
tss_ISelement_a <- left_join(tss_table_antisense,is_pfu_new, by = c("Locus_tag" = "cor_gene"))  %>%
  dplyr::filter(!is.na(id)) %>%
  mutate(aTSS = Pos,
         strand = Strand,
         length = abs(end_read-start_read)) %>%
  dplyr::select(Locus_tag, strand,aTSS, start_read, end_read, length) %>%
  mutate(insertion_strand = ifelse(strand == "+", "-", "+")) %>%
  rowwise() %>%
  mutate(norm = ifelse(insertion_strand == "+" & aTSS < start_read & aTSS > (start_read - 100),
                       0 + ((start_read - aTSS) * (100/100)),
                       ifelse(insertion_strand == "+" & aTSS >= start_read & aTSS < end_read,
                              100 + ((aTSS - start_read) * (100/length)),
                              ifelse(insertion_strand == "+" & aTSS >= end_read & aTSS <= (end_read+100),
                                     200 + ((aTSS - end_read) * (100/100)),
                                     ifelse(insertion_strand == "-" & aTSS > end_read & aTSS < (end_read +  100),
                                            0 + (aTSS - end_read) * (100/100),
                                            ifelse(insertion_strand == "-" & aTSS <= end_read & aTSS > start_read, 
                                                   100 + (end_read - aTSS) * (100/length),
                                                   ifelse(insertion_strand == "-" & aTSS < start_read & aTSS >= (start_read-100), 
                                                          200 + (start_read - aTSS) * (100/100), NA)))))))

tss_ISelement_a2 <- tss_ISelement_a[!duplicated(tss_ISelement_a$aTSS),] %>%
  dplyr::filter(!is.na(norm))
tss_ISelement_a2

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
pdf(here("figures/rnaseq_figures/iselement_antisense.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot(data = tss_ISelement_a2, aes(x = norm, y = ..count..)) +
  theme_Publication() +
  geom_histogram(fill = viridis_pal()(10)[4],      # Histogram with density instead of count on y-axis
                 binwidth=6, size = 0.3, alpha = 1,colour = "white") +
  guides(fill = F) +
  ylab("Number of aTSS") +
  scale_x_continuous(limits = c(0,300), breaks = c(0,100,200,300), labels = c("-100", "start","stop", "+100")) +
  scale_fill_viridis(alpha = 0.3,discrete = F, option = "viridis", begin = 0.2, end = 0.6) +
  xlab("relative position to IS elements") 
dev.off()





