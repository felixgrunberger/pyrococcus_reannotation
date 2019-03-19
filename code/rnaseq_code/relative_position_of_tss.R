# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Position of start site relative to CDS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(here)
library(tidyverse)
library(data.table)
library(viridis)
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
# load data sets
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> TSSpredator mastertable | split according to TSS
mastertable <- here("data/annogesic_data/MasterTable.tsv")

tss_table_primary <- fread(mastertable) %>%
  as.data.table() %>%
  dplyr::filter(Primary == 1) %>%
  mutate(type = "Primary") %>%
  dplyr::select(Locus_tag, Strand, Pos,type, UTRlength) 

tss_table_secondary <- fread(mastertable) %>%
  as.data.table() %>%
  dplyr::filter(Secondary == 1) %>%
  mutate(type = "Secondary") %>%
  dplyr::select(Locus_tag, Strand, Pos,type, UTRlength) 

tss_table_internal <- fread(mastertable) %>%
  as.data.table() %>%
  dplyr::filter(Internal == 1) %>%
  mutate(type = "Internal") %>%
  dplyr::select(Locus_tag, Strand, Pos,type, UTRlength) 

tss_table_antisense <- fread(mastertable) %>%
  as.data.table() %>%
  dplyr::filter(Antisense == 1) %>%
  mutate(type = "Antisense") %>%
  dplyr::select(Locus_tag, Strand, Pos,type, UTRlength) 

# >> load genome file .gff
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

# >> calculate relative position for internal TSS (in a window -300 upstream of start, +300 downstream of stop, and scaled according to gene length in between)
tss_gff_i <- left_join(tss_table_internal,gff_table, by = c("Locus_tag" = "id")) %>%
  dplyr::select(Pos,Locus_tag, strand, start, end, length) %>%
  mutate(norm = ifelse(strand == "+" & Pos <= start & Pos > (start - 300),
                       0 + (Pos - start + 300) * (300/300),
                       ifelse(strand == "+" & Pos > start & Pos <= end,
                              300 + (Pos - start) * (300/length),
                              ifelse(strand == "+" & Pos > end & Pos <= (end+300),
                                     600 + (Pos - end) * (300/300),
                                     ifelse(strand == "-" & Pos >= end & Pos < (end + 300),
                                            0 + (end - Pos + 300) * (300/300),
                                            ifelse(strand == "-" & Pos < end & Pos >= start, 
                                                   300 + (end - Pos) * (300/length),
                                                   ifelse(strand == "-" & Pos < start & Pos >= (start-300), 
                                                          600 + (start - Pos) * (300/300), NA)))))))

# >> calculate relative position for antisense TSS (in a window -300 upstream of start, +300 downstream of stop, and scaled according to gene length in between)
tss_gff_a <- left_join(tss_table_antisense,gff_table, by = c("Locus_tag" = "id")) %>%
  dplyr::select(Pos,Locus_tag, strand, start, end, length) %>%
  mutate(norm = ifelse(strand == "+" & Pos <= start & Pos > (start - 300),
                       0 + (Pos - start + 300) * (300/300),
                       ifelse(strand == "+" & Pos > start & Pos <= end,
                              300 + (Pos - start) * (300/length),
                              ifelse(strand == "+" & Pos > end & Pos <= (end+300),
                                     600 + (Pos - end) * (300/300),
                                     ifelse(strand == "-" & Pos >= end & Pos < (end + 300),
                                            0 + (end - Pos + 300) * (300/300),
                                            ifelse(strand == "-" & Pos < end & Pos >= start, 
                                                   300 + (end - Pos) * (300/length),
                                                   ifelse(strand == "-" & Pos < start & Pos >= (start-300), 
                                                          600 + (start - Pos) * (300/300), NA)))))))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot density curves
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> primary
pdf(here("figures/rnaseq_figures/position_primaryTSS.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot(data = tss_table_primary, aes(x = 300 - UTRlength)) +
  theme_Publication() +
  guides(fill = F) +
  scale_x_continuous(limits = c(0,900), breaks = c(0,300,600,900), labels = c("-300", "start", "stop","+300")) +
  scale_fill_viridis(discrete = F, option = "viridis", begin = 0, end = 0.8) +
  theme(text = element_text(size = 26)) +
  xlab("UTR length") +
  stat_density(adjust = 0.4,alpha=0.6, size = 1.6, 
               fill = viridis_pal(alpha = 0.5,option = "viridis")(50)[15],
               col = viridis_pal(option = "viridis")(50)[15]) +
  xlab("relative position")
dev.off()

# >> secondary
pdf(here("figures/rnaseq_figures/position_secondaryTSS.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot(data = tss_table_secondary, aes(x = 300 - UTRlength)) +
  theme_Publication() +
  guides(fill = F) +
  scale_x_continuous(limits = c(0,900), breaks = c(0,300,600,900), labels = c("-300", "start", "stop","+300")) +
  scale_fill_viridis(discrete = F, option = "viridis", begin = 0, end = 0.8) +
  theme(text = element_text(size = 26)) +
  stat_density(adjust = 0.4,alpha=0.6, size = 1.6, 
               fill = viridis_pal(alpha = 0.5,option = "viridis")(50)[15],
               col = viridis_pal(option = "viridis")(50)[15]) +
  xlab("relative position")
dev.off()

# >> internal
pdf(here("figures/rnaseq_figures/position_internalTSS.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot(data = tss_gff_i, aes(x = norm)) +
  theme_Publication() +
  guides(fill = F) +
  scale_x_continuous(limits = c(0,900), breaks = c(0,300,600,900), labels = c("-300", "start", "stop","+300")) +
  scale_fill_viridis(discrete = F, option = "viridis", begin = 0, end = 0.8) +
  theme(text = element_text(size = 26)) +
  stat_density(adjust = 0.4,alpha=0.6, size = 1.6, 
               fill = viridis_pal(alpha = 0.5,option = "viridis")(50)[15],
               col = viridis_pal(option = "viridis")(50)[15]) +
  xlab("relative position") 
dev.off()

# >> antisense
pdf(here("figures/rnaseq_figures/position_antisenseTSS.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot(data = tss_gff_a, aes(x = norm)) +
  theme_Publication() +
  guides(fill = F) +
  scale_x_continuous(limits = c(0,900), breaks = c(0,300,600,900), labels = c("-300", "start", "stop","+300")) +
  scale_fill_viridis(discrete = F, option = "viridis", begin = 0, end = 0.8) +
  theme(text = element_text(size = 26)) +
  stat_density(adjust = 0.4,alpha=0.6, size = 1.6, 
               fill = viridis_pal(alpha = 0.5,option = "viridis")(50)[15],
               col = viridis_pal(option = "viridis")(50)[15]) +
  xlab("relative position")
dev.off()

