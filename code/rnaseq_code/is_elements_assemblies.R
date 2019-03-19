# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Number of iselements in pyrococcus assemblies
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(here)
library(data.table)
library(viridis)
library(ape)
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
# load ISEScan 1.6 data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> COM1
com1_ises <- read.gff(here("data/ises_data/pfu_com1_ises.fasta.gff")) %>%
  dplyr::filter(type == "insertion_sequence") %>%
  summarise(counts = n(),
            group = "COM1")

# >> DSM3638
dsmz_old <- read.gff(here("data/ises_data/pfu_dsm3638_ises.fasta.gff")) %>%
  dplyr::filter(type == "insertion_sequence") %>%
  summarise(counts = n(),
            group = "DSM 3638 \n(2001)")

# >> DSMZ NEW ASSEMBLY
dsmz_new_assembly <- read.gff(here("data/ises_data/pfu_dsmz_assembly_ises.fasta.gff")) %>%
  dplyr::filter(type == "insertion_sequence") %>%
  summarise(counts = n(),
            group = "DSM 3638 \n(2016)")

# >> combine groups
ises_table <- rbind(com1_ises, dsmz_old, dsmz_new_assembly)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
pdf(here("figures/rnaseq_figures/iselements_pyrococcus_strains.pdf"), 
    width = 5, height = 3.5, paper = "special",onefile=FALSE)

ggplot(data = ises_table, aes(y = counts, x = group, fill = group)) +
  geom_bar(stat="identity") +
  theme_Publication() +
  guides(fill = F) +
  scale_fill_viridis(alpha = 1,discrete = T, option = "viridis", begin = 0.9, end = 0.4) +
  geom_text(aes(label=counts), vjust=1.6, color="white", size=3.5)+
  ylab("number of IS elements") +
  xlab("assembly/strain")

dev.off()




