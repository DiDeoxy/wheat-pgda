library(SNPRelate)
library(plyr)

## load data into R object from the GDS object using the script at the source location
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys.gds"
source("src\\functions\\data_loading.R")

## open the GDS object again and select SNPs with a missing data rate below 0.1
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys.gds")
snpset.ids.list <- snpgdsSelectSNP(wheat, autosome.only = F, missing.rate = 0.1)
snpgdsClose(wheat)
snp.indices <- match(unlist(snpset.ids.list), snp.id)

## create a new gds object with only the snps with less than 10% missing data
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_phys_snpgdsSelectSNP.gds",
                 genmat = genotypes[snp.indices,],
                 sample.id = sample.id, 
                 snp.id = snp.id[snp.indices],
                 snp.chromosome = snp.chrom[snp.indices],
                 snp.position = snp.pos[snp.indices],
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)

## same aas above but for gen map gds
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_gen.gds"
source("src\\functions\\data_loading.R")
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_gen_snpgdsSelectSNP.gds",
                 genmat = genotypes[snp.indices,],
                 sample.id = sample.id, 
                 snp.id = snp.id[snp.indices],
                 snp.chromosome = snp.chrom[snp.indices],
                 snp.position = snp.pos[snp.indices],
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)

## eliminate those individuals that show identity by state (IBS, fractional identity) greater than 0.99
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_snpgdsSelectSNP.gds"
source("src\\functions\\data_loading.R")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_snpgdsSelectSNP.gds")
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

## i take the pairs and use a perl function to find the connected graphs
## i.e. those sets of pairs that form a connected triangle, square or more

## eliminate individuals so that only one from each pair, triangle, or square remains
NILs <- c("PT434", "BW811", "AC Minto 1", "Avocet 1", "BW275 1", 
          "PT616", "BW427 1", "BW492", "BW948", "Carberry 1", 
          "CDC Stanley 1", "PT754", "SWS349", "Somerset 1", "Stettler 1",
          "AC Reed 1", "SWS87", "SWS241", "SWS345", "SWS363",
          "SWS390", "SWS408", "SWS410")

## find the indices of the NILs
sample.indices <- match(NILs, sample.id)

## eliminate the nils from the sample annotation info
retainedSampAnnot <-  list()
for(name in names(samp.annot)) {
  retainedSampAnnot[[name]] <- factor(samp.annot[[name]][-sample.indices])
}
##create gds ubjects without the NILs, physical map
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds",
                 genmat = genotypes[,-sample.indices],
                 sample.id = sample.id[-sample.indices], 
                 snp.id = snp.id,
                 snp.chromosome = snp.chrom,
                 snp.position = snp.pos,
                 other.vars = list(samp.annot = retainedSampAnnot),
                 snpfirstdim = T)

## gen map
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_gen_snpgdsSelectSNP.gds"
source("src\\functions\\data_loading.R")
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_gen_subset_sample.gds",
                 genmat = genotypes[,-sample.indices],
                 sample.id = sample.id[-sample.indices], 
                 snp.id = snp.id,
                 snp.chromosome = snp.chrom,
                 snp.position = snp.pos,
                 other.vars = list(samp.annot = retainedSampAnnot),
                 snpfirstdim = T)

## identify and remove from the overall dataset LD pruned markers
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds"
source("src\\functions\\data_loading.R")

## performs LD pruning to produce a set of SNPs that maximally represent the diversity
## of the genome with as little redundant info as possible used for estimation of relationships
## genome wide (i.e. clustering)
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")
set.seed(1000)
snpset.ids.list <- snpgdsLDpruning(wheat, autosome.only = F, missing.rate = 0.1, maf = 0.05,
                                   ld.threshold = 0.6, slide.max.bp = 1e7)
snpgdsClose(wheat)
snp.indices <- match(unlist(snpset.ids.list), snp.id)

## create a subset gds object containg only the SNPs remaing after ld pruning
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds",
                 genmat = genotypes[snp.indices,],
                 sample.id = sample.id, 
                 snp.id = snp.id[snp.indices],
                 snp.chromosome = snp.chrom[snp.indices],
                 snp.position = snp.pos[snp.indices],
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)

## creating gds subsets for major phenotypic types in the data set
## used in analyses of diversity within these groups
HRS <- which(retainedSampAnnot$designation == "HRS")

load("Data\\Intermediate\\dbscan\\dbscan.RData")
dbscan$cluster[which(dbscan$cluster == 0)] <- dbscan$cluster[which(dbscan$cluster == 0)] + 6
dbscan6 <- factor(dbscan$cluster)

indexNotCanadian <- c(which(dbscan6 != 1))
HRS <- HRS[-which(HRS %in% indexNotCanadian) ]

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