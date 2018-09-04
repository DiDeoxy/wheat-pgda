library(poppr)
library(SNPRelate)
library(ape)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA")

## loading the gds of the data and pullling some attributes out
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds"
source("Analysis\\R\\functions\\data_loading.R")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
distances <- 1 - snpgdsIBS(wheat, autosome.only = F)$ibs
snpgdsClose(wheat)

for (group in c("HRS_SWS", "HRS_HRW", "HRW_SWS")) {
  genind <- get(load(paste0("Data\\Intermediate\\Adegenet\\genind_", group, ".RData")))
  
  genindSeploc <- seploc(genind)
  for (i in 1:length(genindSeploc)){
    strata(genindSeploc[[i]]) <- strata(genind)
  }
  wheatAmova <- list()
  for (i in 1:length(genindSeploc)) {
    if (isPoly(genindSeploc[[i]])) {
      locus <- missingno(genindSeploc[[i]], type = "geno")
      indivs <- match(rownames(locus$tab), sample.id)
      wheatAmova[[i]] <- poppr.amova(locus, hier = as.formula(paste0("~", group)), missing = "genotype", within = F, clonecorrect = F, 
                                     dist = as.dist(distances[indivs, indivs]))
    } else {
      wheatAmova[[i]] <- NA
    }
  }
  
  save(wheatAmova, file = paste0("Data\\Intermediate\\Adegenet\\amova_", group, ".RData"))
}