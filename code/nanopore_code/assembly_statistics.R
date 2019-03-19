# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# calculate statistics from fasta files
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(here)
library(Biostrings)
library(tidyverse)
library(data.table)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load fasta files of different assemblies with Biostrings
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
inFile_fasta    <- readDNAStringSet(here("data/fasta_data/pfu_dsmz_assembly.fasta"))
inFile_polished <- readDNAStringSet(here("data/fasta_data/pfu_nanopolish.fasta"))
inFile_draft    <- readDNAStringSet(here("data/fasta_data/pfu_canu.fasta"))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# calculate length, gc content, contigs, 1-to-1 identity from dnadiff (Mummer output)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> length
dsmz_length <- sum(letterFrequency(inFile_fasta, letters = "ATGC"))
canu_length <- sum(letterFrequency(inFile_draft, letters = "ATGC"))
nanopolish_length <- sum(letterFrequency(inFile_polished, letters = "ATGC"))

# >> gc content
dsmz_gc <- round(sum(letterFrequency(inFile_fasta, letters = "GC"))/dsmz_length*100, digits = 2)
canu_gc <- round(sum(letterFrequency(inFile_draft, letters = "GC"))/canu_length*100, digits = 2)
nanopolish_gc <- round(sum(letterFrequency(inFile_polished, letters = "GC"))/nanopolish_length*100, digits = 2)

# >> contigs
dsmz_contigs <- length(inFile_fasta)
canu_contigs <- length(inFile_draft)
nanopolish_contigs <- length(inFile_polished)

# >> 1-to-1 identity
canu_diff_file <- fread(here("data/dnadiff_data/canu.dnadiff.report"), fill = T) %>%
  as.data.frame() 
nanopolish_diff_file <- fread(here("data/dnadiff_data/nanopolish.dnadiff.report"), fill = T) %>%
  as.data.frame() 

canu_dnadiff <- as.numeric(canu_diff_file[canu_diff_file$V1 == "AvgIdentity",][1,3])
nanopolish_dnadiff <- as.numeric(nanopolish_diff_file[nanopolish_diff_file$V1 == "AvgIdentity",][1,3])


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# save final table as tsv
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
assembly_statistics <- data.frame("assemblies" = c("Illumina/PacBio", "Nanopore raw (canu assembly)", "Nanopore polished (nanopolish assembly)"),
                                  "total bases" = c(dsmz_length, canu_length, nanopolish_length),
                                  "gc content (percentage)" = c(dsmz_gc, canu_gc, nanopolish_gc),
                                  "contigs" = c(dsmz_contigs, canu_contigs, nanopolish_contigs),
                                  "percentage identity (1 to 1 dnadiff)" = c("100", canu_dnadiff, nanopolish_dnadiff))
write.table(assembly_statistics, file = here("figures/nanopore_figures/assembly_table.tsv"), sep = "\t", row.names = F, quote = F)




