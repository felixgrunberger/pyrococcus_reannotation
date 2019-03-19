# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# Visualize number of TSS´s in each category
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
library(tidyverse)
library(here)
library(data.table)
library(viridis)
library(UpSetR)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# load mastertable from TSSpredator
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
mastertable <- fread(here("data/annogesic_data/MasterTable.tsv")) %>%
  as.data.table()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# preparation for plotting 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #

# >> filter for each category and save in new variables
primary <- mastertable %>%
  dplyr::filter(Primary == 1)
secondary <- mastertable %>%
  dplyr::filter(Secondary == 1)
antisense <- mastertable %>%
  dplyr::filter(Antisense == 1)
internal <- mastertable %>%
  dplyr::filter(Internal == 1)

# >> calculate intersects (off all possible variations) and save in new list
upsetr_table <- c(`Primary` = nrow(primary), 
                  `Secondary` = nrow(secondary), 
                  `Antisense` = nrow(antisense), 
                  `Internal` = nrow(internal), 
                  `Primary&Secondary` = length(intersect(primary$SuperPos,secondary$SuperPos)), 
                  `Primary&Antisense` = length(intersect(primary$SuperPos,antisense$SuperPos)), 
                  `Primary&Internal` = length(intersect(primary$SuperPos,internal$SuperPos)), 
                  `Secondary&Antisense` = length(intersect(secondary$SuperPos,antisense$SuperPos)),
                  `Secondary&Internal` = length(intersect(secondary$SuperPos,internal$SuperPos)),
                  `Antisense&Internal` = length(intersect(antisense$SuperPos,internal$SuperPos)),
                  `Primary&Secondary&Antisense` = length(intersect(intersect(primary$SuperPos,secondary$SuperPos),antisense$SuperPos)),
                  `Secondary&Antisense&Internal`=length(intersect(intersect(secondary$SuperPos,antisense$SuperPos),internal$SuperPos)),
                  `Primary&Secondary&Internal`=length(intersect(intersect(primary$SuperPos,secondary$SuperPos),internal$SuperPos)),
                  `Primary&Antisense&Internal`=length(intersect(intersect(primary$SuperPos,antisense$SuperPos),internal$SuperPos)),
                  `Primary&Secondary&Antisense&Internal` = length(intersect(intersect(intersect(primary$SuperPos,secondary$SuperPos),internal$SuperPos),antisense$SuperPos)))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
# plot using UpSetR package
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
pdf(here("figures/rnaseq_figures/tss_categories.pdf"), 
    width = 10, height = 7, paper = "special",onefile=FALSE)

upset(fromExpression(upsetr_table), scale.sets = "identity", 
      sets = rev(c("Primary", "Secondary","Internal", "Antisense")),
      keep.order = TRUE, nsets = 4, line.size = 3, point.size = 11, 
      sets.x.label = "Number of TSS in category", mainbar.y.label = "Intersections", 
      text.scale = c(1.6, 1.2, 1.2, 1.2, 2, 1.6), 
      main.bar.color = viridis_pal(option = "viridis", begin = 0.2, end = 0.6)(9)[1:9], 
      sets.bar.color = "gray60",
      matrix.color = viridis_pal(option = "viridis",begin = 0.2, end = 0.6)(9)[5])

dev.off()