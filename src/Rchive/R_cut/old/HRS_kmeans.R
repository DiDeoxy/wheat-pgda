## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(plyr)
library(ConsensusClusterPlus)
# 
# ## loading the gds of the data and pullling some attributes out
# wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\hrs_phys_subset_both.gds")
# genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
# genotypes <- replace(genotypes, genotypes == 3, NA)
# dist <- as.dist(1 - snpgdsIBS(wheat, autosome.only = F)$ibs)
# snpgdsClose(wheat)
# dim(genotypes)


# dir <- "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\CCP\\hrs"
# results <- ConsensusClusterPlus(dist, maxK = 20, reps = 1000, pItem = 0.8,
#                                 title = dir, clusterAlg = "kmdist", seed = 100, plot = "png")
# 
# save(results, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\hrs_CCP.RDATA")
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\hrs_CCP.RDATA")

icl = calcICL(results, title=dir, plot="png")
clusters <- icl[["itemConsensus"]][which(icl[["itemConsensus"]][,1] == 4),]

bests <- vector(length=234)
for (i in 1:length(bests)) {
  consensus <- clusters[which(clusters[,3] == i),]
  best <- consensus[which.max(consensus[,4]),][,2]
  bests[i] <- best
}

save(bests, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\hrs_kmeans.RDATA")
