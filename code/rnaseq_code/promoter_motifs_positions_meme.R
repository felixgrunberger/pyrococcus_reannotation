# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Promoter motifs (consensus of pTSS & best match of all 4 TSS categories)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(here)
library(data.table)
library(viridis)
library(Biostrings)
library(ggseqlogo)
library(seqinr)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load mastertable from TSSpredator & fasta file
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> fasta
pfu_fasta <- readDNAStringSet(file = here("data/fasta_data/pfu_dsmz_assembly.fasta"))
names(pfu_fasta) <-  "CP023154"

# >> tss table
mastertable <- fread(here("data/annogesic_data/MasterTable.tsv")) %>%
  as.data.table() 


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# calculate PWM for all categories and plot motif
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> strand specific sequences from -50 -> +10
mastertable_all <- mastertable %>%
  rowwise() %>% 
  mutate(cons_seq = ifelse(Strand == "+", as.character(pfu_fasta$CP023154[(Pos-50):(Pos+10)]), 
                           as.character(reverseComplement(pfu_fasta$CP023154[(Pos-10):(Pos+50)])))) %>%
  as.data.table() %>%
  select(cons_seq) 

# >> colors for bases
color_scale = make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                              cols=viridis_pal(begin = 0.2, end = 0.6,option = "viridis")(25)[c(1,7,13,19)])


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot with ggseqlogo
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
pdf(here("figures/rnaseq_figures/pwm_all_motif.pdf"), 
    width = 26, height = 3, paper = "special",onefile=FALSE)
ggplot() + 
  geom_logo(mastertable_all$cons_seq, font = "helvetica_regular", col_scheme = color_scale, seq_type = "dna") + 
  theme_logo() 
dev.off()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# extract sequences and move to meme search for best fitting motif
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> get 50 bases upstream of TSS (strand specific)
tss <- mastertable %>%
  rowwise() %>% 
  mutate(promoter_seq = ifelse(Strand == "+", as.character(pfu_fasta$CP023154[(Pos-50):(Pos)]), 
                               as.character(reverseComplement(pfu_fasta$CP023154[(Pos):(Pos+50)])))) %>%
  as.data.table()

# >> write sequences to fasta file (package seqinr)
write.fasta(as.list(tss$promoter_seq[tss$Primary == 1]), 
            as.list(tss$Pos[tss$Primary == 1]), file = here("data/meme_data/primary_meme.fasta")) 
write.fasta(as.list(tss$promoter_seq[tss$Secondary == 1]), 
            as.list(tss$Pos[tss$Secondary == 1]), file = here("data/meme_data/secondary_meme.fasta")) 
write.fasta(as.list(tss$promoter_seq[tss$Internal == 1]), 
            as.list(tss$Pos[tss$Internal == 1]), file = here("data/meme_data/internal_meme.fasta")) 
write.fasta(as.list(tss$promoter_seq[tss$Antisense == 1]), 
            as.list(tss$Pos[tss$Antisense == 1]), file = here("data/meme_data/antisense_meme.fasta")) 

# >> analyze with meme commandline version (default options except `search given strand only`)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load output from meme
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
primary_motif <- read.table(here("data/meme_data/primary_motif_1_fasta.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1)) 
primary_motif <-   c(as.character(droplevels(primary_motif$V1)))
secondary_motif <- read.table(here("data/meme_data/secondary_motif_1_fasta.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1))
secondary_motif <- c(as.character(droplevels(secondary_motif$V1)))
internal_motif <- read.table(here("data/meme_data/internal_motif_1_fasta.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1))
internal_motif <- c(as.character(droplevels(internal_motif$V1)))
antisense_motif <- read.table(here("data/meme_data/antisense_motif_1_fasta.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1))
antisense_motif <- c(as.character(droplevels(antisense_motif$V1)))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot with ggseqlogo
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> primary
pdf(here("figures/rnaseq_figures/primary_motif.pdf"), 
    width = nchar(primary_motif)[1]/1.3, height = 3, paper = "special",onefile=FALSE)
ggplot() + 
  geom_logo(primary_motif, font = "helvetica_regular", col_scheme = color_scale, seq_type = "dna") + 
  theme_logo() 
dev.off()
# >> secondary
pdf(here("figures/rnaseq_figures/secondary_motif.pdf"), 
    width = nchar(secondary_motif)[1]/1.3, height = 3, paper = "special",onefile=FALSE)
ggplot() + 
  geom_logo(secondary_motif, font = "helvetica_regular", col_scheme = color_scale, seq_type = "dna") + 
  theme_logo() 
dev.off()
# >> internal
pdf(here("figures/rnaseq_figures/internal_motif.pdf"), 
    width = nchar(internal_motif)[1]/1.3, height = 3, paper = "special",onefile=FALSE)
ggplot() + 
  geom_logo(internal_motif, font = "helvetica_regular", col_scheme = color_scale, seq_type = "dna") + 
  theme_logo() 
dev.off()
# >> antisense
pdf(here("figures/rnaseq_figures/antisense_motif.pdf"), 
    width = nchar(antisense_motif)[1]/1.3, height = 3, paper = "special",onefile=FALSE)
ggplot() + 
  geom_logo(antisense_motif, font = "helvetica_regular", col_scheme = color_scale, seq_type = "dna") + 
  theme_logo() 
dev.off()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# position of logos
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
primary_pos <- fread(here("data/meme_data/primary_position.tsv"), fill = T) %>%
  mutate(position = as.numeric(V2) - 40) %>%
  select(position) 
secondary_pos <- fread(here("data/meme_data/secondary_positions.tsv"), fill = T) %>%
  mutate(position = as.numeric(V2) - 40) %>%
  select(position) 
internal_pos <- fread(here("data/meme_data/internal_positions.tsv"), fill = T) %>%
  mutate(position = as.numeric(V2) - 40) %>%
  select(position) 
antisense_pos <- fread(here("data/meme_data/antisense_positions.tsv"), fill = T) %>%
  mutate(position = as.numeric(V2) - 40) %>%
  select(position) 


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot position of logos aligned to third T of TATA box
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> primary motif position
pdf(here("figures/rnaseq_figures/promoters_primary_position.pdf"), 
    width = 24, height = 2.5, paper = "special",onefile=FALSE)
ggplot(data = primary_pos, aes(x = position, y = ..count..)) +
  theme_light() +
  geom_histogram(aes(y=..count.., fill = ..count..),      # Histogram with density instead of count on y-axis
                 binwidth=1, size = 0.3, alpha = 0.6,colour = "white") +
  scale_fill_viridis(alpha = 0.3,discrete = F, option = "viridis", begin = 0.2, end = 0.6) +
  guides(fill = F) +
  ylab("number") +
  xlab("") +
  theme(text = element_text(size = 26))
dev.off()  


# >> secondary motif position 
pdf(here("figures/rnaseq_figures/promoters_secondary_position.pdf"), 
    width = 24, height = 2.5, paper = "special",onefile=FALSE)
ggplot(data = secondary_pos, aes(x = position, y = ..count..)) +
  theme_light() +
  geom_histogram(aes(y=..count.., fill = ..count..),      # Histogram with density instead of count on y-axis
                 binwidth=1, size = 0.3, alpha = 0.6,colour = "white") +
  scale_fill_viridis(alpha = 0.3,discrete = F, option = "viridis", begin = 0.2, end = 0.6) +
  guides(fill = F) +
  ylab("number") +
  xlab("") +
  theme(text = element_text(size = 26))
dev.off()  


# >> internal motif position 
pdf(here("figures/rnaseq_figures/promoters_internal_position.pdf"), 
    width = 24, height = 2.5, paper = "special",onefile=FALSE)
ggplot(data = internal_pos, aes(x = position, y = ..count..)) +
  theme_light() +
  geom_histogram(aes(y=..count.., fill = ..count..),      # Histogram with density instead of count on y-axis
                 binwidth=1, size = 0.3, alpha = 0.6,colour = "white") +
  scale_fill_viridis(alpha = 0.3,discrete = F, option = "viridis", begin = 0.2, end = 0.6) +
  guides(fill = F) +
  ylab("number") +
  xlab("") +
  theme(text = element_text(size = 26))
dev.off()  


# >> antisense motif position 
pdf(here("figures/rnaseq_figures/promoters_antisense_position.pdf"), 
    width = 24, height = 2.5, paper = "special",onefile=FALSE)
ggplot(data = antisense_pos, aes(x = position, y = ..count..)) +
  theme_light() +
  geom_histogram(aes(y=..count.., fill = ..count..),      # Histogram with density instead of count on y-axis
                 binwidth=1, size = 0.3, alpha = 0.6,colour = "white") +
  scale_fill_viridis(alpha = 0.3,discrete = F, option = "viridis", begin = 0.2, end = 0.6) +
  guides(fill = F) +
  ylab("number") +
  xlab("") +
  theme(text = element_text(size = 26))
dev.off()  