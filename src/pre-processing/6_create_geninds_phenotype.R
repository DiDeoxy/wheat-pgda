library(tidyverse)
library(SNPRelate)
library(adegenet)

## loading the gds of the data and pullling some attributes out
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_both.gds"
source("src\\functions\\data_loading.R")

# making the data palatable by genind
genotypes <- geontypes %>%
             replace(. == 0, "A") %>%
             replace(. == 2, "B") %>%
             replace(. == 3, "")

# for making subsets of the samples for the different phenotypes
index_hard <- c(which(desig == "HRS"), which(desig == "HWS"), 
               which(desig == "HRW"), which(desig == "HWW")) # HRS HRW
index_soft <- c(which(desig == "SRS"), which(desig == "SWS"),
               which(desig == "SRW"), which(desig == "SWW")) # SWS
index_red <- c(which(desig == "HRS"), which(desig == "SRS"), 
              which(desig == "HRW"), which(desig == "SRW")) # HRS HRW
index_white <- c(which(desig == "HWS"), which(desig == "SWS"), 
                which(desig == "HWW"), which(desig == "SWW")) # SWS
index_spring <- c(which(desig == "HRS"), which(desig == "HWS"), 
                 which(desig == "SRS"), which(desig == "SWS")) # HRS SWS
index_winter <- c(which(desig == "HRW"), which(desig == "HWW"),
                 which(desig == "SRW"), which(desig == "SWW")) # HRW

## make unknown groups uniformaly NA
era[which(era == "UNKNOWN")] <- "NA"
bp[c(which(bp == "FOREIGN"), which(bp == "N/A"), which(bp == "USA"))] <- "NA"
desig[which(desig == "UNKNOWN")] <- "NA"

# create a hard/soft factor
desig_hs <- as.character(desig)
desig_hs[index_hard] <- "hard"
desig_hs[index_soft] <- "soft"
desig_hs[-c(index_hard, index_soft)] <- "NA"

# create a red/white factor
desig_rw <- as.character(desig)
desig_rw[index_red] <- "red"
desig_rw[index_white] <- "white"
desig_rw[-c(index_red, index_white)] <- "NA"

# create a spring/winter factor
desig_sw <- as.character(desig)
desig_sw[index_winter] <- "winter"
desig_sw[index_spring] <- "spring"
desig_sw[-c(index_winter, index_spring)] <- "NA"

# create combinations of hard/soft, red/white, and spring/winter
desig_sw_rw <- paste(desig_sw, desig_rw)
desig_sw_hs <- paste(desig_sw, desig_hs)
desig_rw_hs <- paste(desig_sw, desig_rw)

# create a hrs vs all other group
desig_hrs_other <- as.character(desig)
desig_hrs_other[which(desig_hrs_other != "HRS")] <- "NA"

# create a data frame of the different strata
strata <- data.frame("SW" = as.factor(desig_sw), "RW" = as.factor(desig_rw),
                     "HS" = as.factor(desig_hs), 
                     "SWRW" = as.factor(desig_sw_rw),
                     "SWHS" = as.factor(desig_sw_hs),
                     "RWHS" = as.factor(desig_rw_hs),
                     "HRSOther" = as.factor(desig_hrs_other),
                     "desig" = as.factor(desig), "era" = as.factor(era),
                     "bp"= as.factor(bp))
                    #  , "hdbscan2" = wheatHdbscan2, 
                    #  "hdbscan6" = wheatHdbscan6, "dend" = clusters)

# turn the genotype data and strate into a genind
genind <- df2genind(t(data.frame(genotypes)),
                    ind.names = as.character(sample_id), loc.names = snp_id,
                    ploidy = 1, type = "codom", ncode = 1, strata = strata)

save(genind, file = paste0("Data\\Intermediate\\Adegenet\\genind_wheat.RData"))