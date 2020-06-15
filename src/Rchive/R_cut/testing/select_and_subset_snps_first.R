

suppressMessages(library(SNPRelate))
library(plyr)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## eliminate those SNPs with a missing rate above 0.10
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys.gds"
source("Analysis\\R\\functions\\data_loading.R")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys.gds")
snpset.ids.list <- snpgdsSelectSNP(wheat, autosome.only = F, missing.rate = 0.1)
snpgdsClose(wheat)
snp.indices <- match(unlist(snpset.ids.list), snp.id)

snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_phys_snpgdsSelectSNP.gds",
                 genmat = genotypes[snp.indices,],
                 sample.id = sample.id, 
                 snp.id = snp.id[snp.indices],
                 snp.chromosome = snp.chrom[snp.indices],
                 snp.position = snp.pos[snp.indices],
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)

## eliminate those individuals that show identity by state (IBS, fractional identity) greater than 0.995
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_snpgdsSelectSNP.gds"
source("Analysis\\R\\functions\\data_loading.R")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_snpgdsSelectSNP.gds")
set.seed(1000)
snpset.ids.list <- snpgdsLDpruning(wheat, autosome.only = F, ld.threshold = 0.6, 
                                   slide.max.bp = 1e7)
snpgdsClose(wheat)
snp.indices <- match(unlist(snpset.ids.list), snp.id)

snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_phys_subset_snp.gds",
                 genmat = genotypes[snp.indices,],
                 sample.id = sample.id, 
                 snp.id = snp.id[snp.indices],
                 snp.chromosome = snp.chrom[snp.indices],
                 snp.position = snp.pos[snp.indices],
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)

## identify highly similar samples
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_snp.gds")
IBS <- snpgdsIBS(wheat, autosome.only = F)
snpgdsClose(wheat)

pairs <- which(IBS$ibs >= 0.99, arr.ind = T)

indices <- vector()
for (i in 1:dim(pairs)[1]) {
  if (pairs[i,1] == pairs[i,2]){
    indices <<- c(indices, i)
  }
}
pairs <- pairs[-indices,]
pairs

# NILs <- c("BW811", "AC Minto 1", "AC Reed 1", "Avocet 1", "BW275 1", 
#           "BW427 1", "BW492", "Carberry 1", "CDC Stanley 1", "Somerset 1",
#           "Stettler 1", "SWS408", "SWS241", "SWS345", "SWS410")

NILs <- c("PT434", "BW811", "AC Minto 1", "Avocet 1", "BW275 1", 
          "BW395", "BW427 1", "BW492", "BW922", "BW948",
          "Carberry 1", "CDC Stanley 1", "PT754", "SWS349", "Somerset 1", 
          "Stettler 1", "SWS241", "SWS345", "AC Reed 1", "SWS87", 
          "SWS390", "SWS408", "SWS410")

## create subset of both samples and snps
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_snp.gds"
source("Analysis\\R\\functions\\data_loading.R")

sample.indices <- match(NILs, sample.id)

retainedSampAnnot <-  list()
for(name in names(samp.annot)) {
  retainedSampAnnot[[name]] <- factor(samp.annot[[name]][-sample.indices])
}
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds",
                 genmat = genotypes[,-sample.indices],
                 sample.id = sample.id[-sample.indices], 
                 snp.id = snp.id,
                 snp.chromosome = snp.chrom,
                 snp.position = snp.pos,
                 other.vars = list(samp.annot = retainedSampAnnot),
                 snpfirstdim = T)

## create subsets of samples
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_snpgdsSelectSNP.gds"
source("Analysis\\R\\functions\\data_loading.R")

sample.indices <- match(NILs, sample.id)

retainedSampAnnot <-  list()
for(name in names(samp.annot)) {
  retainedSampAnnot[[name]] <- factor(samp.annot[[name]][-sample.indices])
}
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds",
                 genmat = genotypes[,-sample.indices],
                 sample.id = sample.id[-sample.indices], 
                 snp.id = snp.id,
                 snp.chromosome = snp.chrom,
                 snp.position = snp.pos,
                 other.vars = list(samp.annot = retainedSampAnnot),
                 snpfirstdim = T)

## Subsets by gross phenotypic type
HRS <- which(retainedSampAnnot$designation == "HRS")
HRSRetainedSampAnnot <- list()
for(name in names(samp.annot)) {
  HRSRetainedSampAnnot[[name]] <- factor(retainedSampAnnot[[name]][HRS])
}
snpgdsCreateGeno("Data\\Intermediate\\GDS\\HRS_phys_subset_sample.gds",
                 genmat = genotypes[,HRS],
                 sample.id = sample.id[HRS], 
                 snp.id = snp.id,
                 snp.chromosome = snp.chrom,
                 snp.position = snp.pos,
                 other.vars = list(samp.annot = HRSRetainedSampAnnot),
                 snpfirstdim = T)

HRW <- which(retainedSampAnnot$designation == "HRW")
HRWRetainedSampAnnot <- list()
for(name in names(samp.annot)) {
  HRWRetainedSampAnnot[[name]] <- factor(retainedSampAnnot[[name]][HRW])
}
snpgdsCreateGeno("Data\\Intermediate\\GDS\\HRW_phys_subset_sample.gds",
                 genmat = genotypes[,HRW],
                 sample.id = sample.id[HRW], 
                 snp.id = snp.id,
                 snp.chromosome = snp.chrom,
                 snp.position = snp.pos,
                 other.vars = list(samp.annot = HRWRetainedSampAnnot),
                 snpfirstdim = T)

SWS <- which(retainedSampAnnot$designation == "SWS")
SWSRetainedSampAnnot <- list()
for(name in names(samp.annot)) {
  SWSRetainedSampAnnot[[name]] <- factor(retainedSampAnnot[[name]][SWS])
}
snpgdsCreateGeno("Data\\Intermediate\\GDS\\SWS_phys_subset_sample.gds",
                 genmat = genotypes[,SWS],
                 sample.id = sample.id[SWS], 
                 snp.id = snp.id,
                 snp.chromosome = snp.chrom,
                 snp.position = snp.pos,
                 other.vars = list(samp.annot = SWSRetainedSampAnnot),
                 snpfirstdim = T)