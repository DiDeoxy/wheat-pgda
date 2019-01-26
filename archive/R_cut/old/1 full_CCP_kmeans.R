## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(plyr)
library(ConsensusClusterPlus)
library(adegenet)
library(poppr)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the gds of the data and pullling some attributes out
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
dist <- as.dist(1 - snpgdsIBS(wheat, autosome.only = F)$ibs)
snpgdsClose(wheat)

# genind <- get(load("Data\\Intermediate\\Adegenet\\genind_subset_sample_maf05.RData"))
# dist <- diss.dist(genind)

dir <- "Results\\CCP\\"
results <- ConsensusClusterPlus(dist, maxK = 6, reps = 5000, pItem = 0.8,
                                title = dir, clusterAlg = "kmdist", seed = 100, plot = "png")

save(results, file = "Data\\Intermediate\\CCP\\CCP.RDATA")
load("Data\\Intermediate\\CCP\\CCP.RDATA")
# str(results[[6]])

genotypes <- replace(genotypes, genotypes == 3, NA)

icl = calcICL(results, title=dir, plot="png")
clusters <- icl[["itemConsensus"]][which(icl[["itemConsensus"]][,1] == 6),]
bests <- vector(length=ncol(genotypes))
for (i in 1:ncol(genotypes)) {
  consensus <- clusters[which(clusters[,3] == i),]
  bests[i] <- consensus[which.max(consensus[,4]),][,2]
}
save(bests, file = "Data\\Intermediate\\CCP\\full_kmeans.RDATA")
