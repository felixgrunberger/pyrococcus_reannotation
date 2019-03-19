# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# comparison of nucleotide frequencies (6mer) in pyrococcus assemblies 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(here)
library(data.table)
library(Biostrings)
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
# load fasta files of different assemblies
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
inFile_fasta    <- here("data/fasta_data/pfu_dsmz_assembly.fasta")
inFile_polished <- here("data/fasta_data/pfu_nanopolish.fasta")
inFile_draft    <- here("data/fasta_data/pfu_canu.fasta")

# calculate sixmer frequencies
f6_reference <- Biostrings::oligonucleotideFrequency(readDNAStringSet(filepath = inFile_fasta),6)
f6_polished <- Biostrings::oligonucleotideFrequency(readDNAStringSet(filepath = inFile_polished),6)
f6_draft <- Biostrings::oligonucleotideFrequency(readDNAStringSet(filepath = inFile_draft),6)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# calculate data tables and prepare for plotting
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
reference <- f6_reference %>%
  t() %>%
  as.data.table() %>%
  mutate(n_ref = V1,
         sixmer = unlist(colnames(f6_reference))) %>%
  dplyr::select(sixmer,n_ref)
#
draft <- f6_draft %>%
  t() %>%
  as.data.table() %>%
  mutate(n_draft = colSums(f6_draft)) %>%
  dplyr::select(n_draft)
#
polished <- f6_polished %>%
  t() %>%
  as.data.table() %>%
  mutate(n_polished = colSums(f6_polished)) %>%
  dplyr::select(n_polished)
#
joined <- cbind(reference,draft,polished) 
#
for (i in 1: length(joined$sixmer)){
  joined$nr_cons_letter[i]  <- max(with(rle(strsplit(joined$sixmer[i],"")[[1]]), paste(lengths)))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot (and save as PDF)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
pdf(here("figures/nanopore_figures/sixmer_comparison.pdf"), 
    width = 10, height = 10, paper = "special",onefile=FALSE)
ggplot(data = joined, aes(x = n_ref,fill = as.factor(nr_cons_letter), text = sixmer)) +
  geom_point(aes(y = n_draft),size = 18, alpha = 0.3, shape = 21, fill = "grey",col = "white", stroke = .4) +
  geom_point(aes(y = n_polished),size = 18, alpha = 0.8, shape = 21, col = "white", stroke = .4) +
  theme_Publication() +
  scale_fill_viridis_d() +
  xlab("sixmer counts in reference genome (PacBio/Illumina)") +
  ylab("sixmer counts in polished de novo assembly (MinION)") +
  geom_abline(col = "white", lty = "dashed", lwd = 1) +
  coord_fixed(ratio = 1) +
  labs(fill = "Number of \nconsecutive bases")
dev.off()

