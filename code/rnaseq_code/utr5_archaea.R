# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Archaeal data sets, UTR length, dna gc contents
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(here)
library(data.table)
library(pdftools)
library(stringr)
library(readxl)
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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load archaeal data sets 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> METHANOCALDOCOCCUS JANASCHII
### file location "https://media.nature.com/original/nature-assets/nmicrobiol/2017/nmicrobiol201721/extref/nmicrobiol201721-s3.xlsx",
mja_file <- read_excel("..file_location..", 
                       sheet = 1, col_names = T, col_types = NULL, na = "", skip = 30) %>%
  as.data.table() %>%
  mutate(UTR = `Leader length9`,
         ORG = "MJA",
         MED = median(UTR,na.rm = T)) %>%
  dplyr::filter(!is.na(UTR) & UTR >= 0) %>%
  dplyr::select(UTR, ORG, MED) 

# >> THERMOCOCCUS ONNURINEUS
### file location "https://media.nature.com/original/nature-assets/srep/2017/170220/srep43044/extref/srep43044-s1.pdf"
ton_file1 <- pdf_text("https://media.nature.com/original/nature-assets/srep/2017/170220/srep43044/extref/srep43044-s1.pdf") %>%
  strsplit(split = "\n")
ton_file2 <- unlist(ton_file1)
ton_file <- as.data.table(str_split(ton_file2, " +", n = 20, simplify = T)) %>%
  dplyr::filter(substr(V1,1,2) == "TU" & V5 == "P") %>%
  mutate(V13 = ifelse(V13 == "", NA, V13),
         V12 = ifelse(V12 == "", NA, V12),
         V11 = ifelse(V11 == "", NA, V11),
         UTR = as.numeric(ifelse(!is.na(V13), V13, 
                                 ifelse(!is.na(V12), V12,
                                        ifelse(!is.na(V11), V11,V10))))) %>%
  mutate(UTR = ifelse(UTR < 0, 0, UTR),
         ORG = "TON",
         MED = median(UTR,na.rm = T)) %>%
  dplyr::select(UTR, ORG, MED) 

# >> HALOFERAX VOLCANII
### file location "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-016-2920-y/MediaObjects/12864_2016_2920_MOESM2_ESM.xlsx"
hvo_file <- read_excel("..file_location..") %>%
  mutate(UTR = ifelse(!is.na(`leaderless start distance`),`leaderless start distance`,as.numeric(`UTR Length`)),
         ORG = "HVO",
         MED = median(UTR,na.rm = T)) %>%
  dplyr::select(UTR, ORG, MED) 

# >> THERMOCOCCUS KODAKARENSIS
### file location "https://static-content.springer.com/esm/art%3A10.1186%2F1471-2164-15-684/MediaObjects/12864_2014_6679_MOESM2_ESM.xlsx"
tko_file <- read_excel("..file_location..") %>%
  dplyr::filter(`TSS type` == "primary") %>%
  mutate(UTR1 = ifelse(`Gene strand` == "+",`Gene start`-`TSS position`, `TSS position`-`Gene end`),
         UTR = ifelse(UTR1 < 0, 0 , UTR1),
         ORG = "TKO",
         MED = median(UTR,na.rm = T)) %>%
  dplyr::select(UTR, ORG, MED) 

# >> PYROCOCCUS FURIOSUS
pfu_file <- fread(here("data/annogesic_data/MasterTable.tsv")) %>%
  as.data.table() %>%
  dplyr::filter(Primary == 1) %>%
  mutate(UTR = UTRlength,
         ORG = "PFU",
         MED = median(UTR,na.rm = T)) %>%
  dplyr::filter(!is.na(UTR)) %>%
  dplyr::select(UTR, ORG, MED) 

# >> combine tables
joined_table <- rbind(mja_file,tko_file, hvo_file, ton_file, pfu_file) %>%
  mutate(ORG = reorder(ORG, MED))

joined_table$ORG <- factor(joined_table$ORG,levels(joined_table$ORG)[c(1,2,3,5,4)])


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
joined_table2 <- joined_table %>%
  mutate(UTR2 = UTR + 1)

pdf(here("figures/rnaseq_figures/utr5_archaea.pdf"), 
    width = 8, height = 7, paper = "special",onefile=FALSE)
ggplot(data = joined_table2, aes(x = UTR2, fill = ORG, linetype = ORG)) +
  geom_density(size = 1.2, alpha = 0.6, color = "black") + 
  scale_x_continuous(limits = c(1,300),breaks = c(1,6,11,51,101,301),labels = c(0,5,10,50,100,300),trans = "log10") +
  scale_fill_viridis_d(begin = 0.2, end = 0.8) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed", "dotted", "solid")) +
  theme_Publication() +
  xlab("5`UTR length [nt, log10 scale]") +
  ylab("density")
dev.off()


