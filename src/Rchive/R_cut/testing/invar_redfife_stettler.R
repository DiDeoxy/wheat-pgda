suppressMessages(library(SNPRelate))

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds"
source("Analysis\\R\\functions\\data_loading.R")
genotypes <- replace(genotypes, genotypes == 3, NA)

# wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")
# snpgdsClose(wheat)

HRS <- which(samp.annot$designation == "HRS")

load("Data\\Intermediate\\CCP\\full_kmeans.RDATA")
bests <- factor(bests)
levels(bests) <- c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6")
cadHRS <- which(bests == "Group 1" | bests == "Group 2" | bests == "Group 4")

load("Data\\Intermediate\\CCP\\CCP.RDATA")
labelOrder <- results[[6]]$consensusTree$order

# length(which(rowSums(genotypes[,cadHRS], na.rm = T) == 0 | rowSums(genotypes[,cadHRS], na.rm = T) == 2*length(cadHRS)))
# length(which(rowSums(genotypes, na.rm = T) == 0 | rowSums(genotypes, na.rm = T) == 2*ncol(genotypes)))

# which(apply(genotypes, 1, function (x) { var(x, na.rm = T) }) == 0)



# which(sample.id[labelOrder][cadHRS] == "BW389")
# bests[which(sample.id == "Journey")]
# bests[which(sample.id == "BW839")]
# bests[which(sample.id == "5601HR")]

# dim(genotypes)
# which(sample.id == "Stettler")
# which(sample.id == "Red Fife")

invariant <- which(rowSums(genotypes[,c(302, 281)], na.rm = T) == 0 | rowSums(genotypes[,c(302, 281)], na.rm = T) == 4)
length(invariant)

labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
labels <- as.vector(t(outer(labels, c("_Part1", "_Part2"), paste, sep = "")))

invarData <- do.call(rbind, list(labels[snp.chrom], snp.pos, genotypes[, 302], genotypes[, 281]))
colnames(invarData) <- snp.id
rownames(invarData) <- c("snp.chrom", "snp.pos", "Stettler", "Red Fife")
invarData[1:4, 1:10]
write.csv(invarData, file = "Data\\Intermediate\\stettler_redFife.csv")
