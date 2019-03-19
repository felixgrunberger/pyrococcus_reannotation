# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Bidirectional transcripts, intergene distance
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(here)
library(data.table)
library(rtracklayer)
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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load data and filter
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

#  >>filter for primary TSS
tss_primary <- fread(here("data/annogesic_data/MasterTable.tsv")) %>%
  as.data.table() %>%
  filter(Primary == 1) %>%
  select(Pos, Strand, Locus_tag, Primary)

# >> load gff annotation file and add information about head-to-head/non head-to-head orientation
gff_table <- readGFF(file = here("data/genome_data/CP023154.gff")) %>% 
  as.data.frame() %>%
  dplyr::filter(type == "exon") %>%
  mutate(id = substr(Parent, 1,15)) %>%
  select(strand, id,start,end) %>%
  filter(isUnique(id)) %>%
  mutate(orientation = ifelse(strand == "+" & lag(strand) == "+", "head-to-tail",
                              ifelse(strand == "-" & lead(strand) == "-", "head-to-tail", 
                                     ifelse(strand == "+" & lag(strand) == "-", "head-to-head",
                                            ifelse(strand == "-" & lead(strand) == "+", "head-to-head", NA)))),
         lag_start = ifelse(orientation == "head-to-head" & strand == "+" , lag(end),
                            ifelse(orientation == "head-to-head" & strand == "-", lead(start), NA)),
         lag_end = ifelse(strand == "+" & orientation == "head-to-tail", lag(end),
                          ifelse(strand == "-" & orientation == "head-to-tail", lead(start), NA))) %>%
  select(strand, orientation, id, lag_end, lag_start, start, end)

# >> join 2 tables and filter for orientation and strand
tss_gff_distance <- left_join(gff_table,tss_primary, by = c("id" = "Locus_tag")) %>%
  mutate(lead_start = lead(start)) %>%
  mutate(distance = ifelse(orientation == "head-to-head" & strand == "+", start - lag_start,
                           ifelse(orientation == "head-to-head" & strand == "-", lag_start - end, 
                                  ifelse(orientation == "head-to-tail" & strand == "+", start - lag_end,
                                         ifelse(orientation == "head-to-tail" & strand == "-", lead_start - end,NA))))) %>%
  mutate(chr = "CP023154",
         start = start,
         end = end,
         name = id) %>%
  dplyr::select(Primary,chr,start,end,name, distance, orientation, lag_end, lag_start, strand,Pos)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot distances using ggplot2
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
pdf(here("figures/rnaseq_figures/intergene_distance.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)
ggplot(data = tss_gff_distance[!is.na(tss_gff_distance$orientation),], 
       aes(x = distance,fill = orientation, color = orientation)) +
  geom_density(alpha = 0.7) +
  scale_x_continuous(limits = quantile(tss_gff_distance$distance[!is.na(tss_gff_distance$orientation)], c(0, 0.95))) +
  scale_fill_grey(start = 0, end = 0.6) +
  scale_color_grey(start = 0, end = 0.6) +#viridis(discrete = T, begin = 0.5, end = 0.9, option = "inferno") +
  #scale_color_viridis(discrete = T, begin = 0.2, end = 0.6) + 
  xlab("intergene distance") +
  theme_Publication() 
dev.off()
