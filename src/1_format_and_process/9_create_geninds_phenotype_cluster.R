library(tidyverse)
library(SNPRelate)
library(adegenet)
library(plyr)

## loading the gds of the data and pullling some attributes out
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned_floor.gds"
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
hrs <- as.character(desig)
hrs[which(hrs != "HRS")] <- "Other"

sws <- as.character(desig)
sws[which(sws != "SWS")] <- "Other"

hrw <- as.character(desig)
hrw[which(hrw != "HRW")] <- "Other"

# create each cluster vs all other
cluster1 <- cluster
cluster1[which(cluster1 != 1)] <- 0

cluster2 <- cluster
cluster2[which(cluster2 != 2)] <- 0

cluster3 <- cluster
cluster3[which(cluster3 != 3)] <- 0

cluster4 <- cluster
cluster4[which(cluster4 != 4)] <- 0

cluster5 <- cluster
cluster5[which(cluster5 != 5)] <- 0

# create CSWS

# create a data frame of the different strata
strata <- data.frame(era = era, bp = bp,
                     texture = texture, colour = colour, habit = habit,
                     phenotype = desig,
                     hrs = as.factor(hrs),
                     sws = as.factor(sws),
                     hrw = as.factor(hrw),
                     clusters = as.factor(cluster),
                     cluster1 = as.factor(cluster1),
                     cluster2 = as.factor(cluster2),
                     cluster3 = as.factor(cluster3),
                     cluster4 = as.factor(cluster4),
                     cluster5 = as.factor(cluster5))

# turn the genotype data and strate into a genind
full_genind <- df2genind(t(data.frame(genotypes)),
                         ind.names = as.character(sample_id),
                         loc.names = snp_id, ploidy = 1, type = "codom",
                         ncode = 1, strata = strata)

write_rds(full_genind,
          path = paste0("Data\\Intermediate\\Adegenet\\full_genind.rds"))