library(tidyverse)
library(SNPRelate)
library(adegenet)
library(plyr)

## loading the gds of the data and pullling some attributes out
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds"
source("src\\R_functions\\data_loading.R")
cluster <- read_rds("Data\\Intermediate\\dbscan\\full_hdbscan.rds")$cluster

# making the data palatable by genind
genotypes <- genotypes %>%
             replace(. == 0, "A") %>%
             replace(. == 2, "B") %>%
             replace(. == 3, "")

## make unknown and small groups uniformaly NA
bp <- revalue(bp, replace = (c("FOREIGN" = "N/A", "USA" = "N/A")))

# create a hrs vs all other group
hrs_other <- as.character(desig)
hrs_other[which(hrs_other != "HRS")] <- "Other"
hrs_other <- as.factor(hrs_other)

# create chrs vs all others
chrs_other <- cluster
chrs_other[which(chrs_other != 5)] <- "Other"
chrs_other[which(chrs_other == 5)] <- "CHRS"
chrs_other <- as.factor(chrs_other)

# create a data frame of the different strata
strata <- data.frame(desig = desig, era = era, bp = bp,
                     texture = texture, colour = colour, habit = habit,
                     texture_colour = paste(texture, colour),
                     texture_habit = paste(texture, habit),
                     colour_habit = paste(colour, habit),
                     hrs_other = hrs_other, cluster = as.factor(cluster),
                     chrs_other = chrs_other)

# turn the genotype data and strate into a genind
full_genind <- df2genind(t(data.frame(genotypes)),
                         ind.names = as.character(sample_id),
                         loc.names = snp_id, ploidy = 1, type = "codom",
                         ncode = 1, strata = strata)

write_rds(full_genind,
         path = paste0("Data\\Intermediate\\Adegenet\\full_genind.rds"))