library(tidyverse)
library(poppr)
library(SNPRelate)
library(ape)

## loading the gds of the data and pullling some attributes out
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\functions\\data_loading.R")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys_subset_both.gds")
distances <- 1 - snpgdsIBS(wheat, autosome.only = F)$ibs
snpgdsClose(wheat)

for (group in c("CHRS_CSWS", "HRS_HRW", "HRW_CSWS")) {
    genind <- get(
        load(paste0("Data\\Intermediate\\Adegenet\\genind_", group, ".RData")))
  
    genind_seploc <- seploc(genind)
    for (i in 1:length(genind_seploc)) {
        strata(genind_seploc[[i]]) <- strata(genind)
    }
    subset_amovas <- list()
    for (i in 1:length(genind_seploc)) {
        if (isPoly(genind_seploc[[i]])) {
        locus <- missingno(genind_seploc[[i]], type = "geno")
        indivs <- match(rownames(locus$tab), sample_id)
        subset_amovas[[i]] <- poppr.amova(
            locus, hier = as.formula(paste0("~", group)), missing = "genotype",
            within = F, clonecorrect = F, 
            dist = as.dist(distances[indivs, indivs]))
        } else {
        subset_amovas[[i]] <- NA
        }
    }
    
    write_rds(subset_amovas,
         path = paste0("Data\\Intermediate\\Adegenet\\amova_", group, ".rds"))
}