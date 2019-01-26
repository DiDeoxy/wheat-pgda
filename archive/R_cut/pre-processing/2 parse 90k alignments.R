library(ape)
library(plyr)

chrs <- paste0("chr", as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep=""))))

alignments <- data.frame()
for (chr in chrs) {
  chrFeats <- read.gff(paste0("Data\\Raw\\Maps\\90K_RefSeqv1_physical_maps\\Infinium90K-", chr, ".gff3"))
  # print(nrow(chrFeats))
  parsedPos <- t(apply(chrFeats, 1, function (feat) {
    attributes <- unlist(strsplit(feat[9], split = ";"))
    coverage = unlist(strsplit(attributes[1], split = "="))[2]
    name = unlist(strsplit(attributes[2], split = "="))[2]
    identity = unlist(strsplit(attributes[3], split = "="))[2]
    ID = unlist(strsplit(attributes[4], split = "="))[2]
    return(c(name, ID, coverage, identity, substr(feat[1], 4, nchar(feat[1])), floor((as.integer(feat[4]) + as.integer(feat[5]))/2), use.names = F))
  }))
  alignments <- rbind(alignments, parsedPos, stringsAsFactors = F)
}

wangGen <- read.csv("Data\\Intermediate\\Maps\\Phys\\pozniak_filtered_map.csv",
                    header = T, col.names = c("Name", "Contig", "Position"), row.names = 1, stringsAsFactors = F)

physMap <- by(alignments, as.factor(alignments[,1]), function (marker) {
  bestAlignments <- data.frame()
  for (i in 1:nrow(marker)) {
    if (as.integer(marker[i,3]) >= 90 & as.integer(marker[i,4]) >= 98) {
      bestAlignments <- rbind(bestAlignments, marker[i,])
    }
  }
  if (nrow(bestAlignments) == 0 | nrow(bestAlignments) > 3) {
    return()
  } else if (nrow(bestAlignments) == 1 & bestAlignments[1,1] %in% row.names(wangGen) & wangGen[bestAlignments[1,1],1] == bestAlignments[1,5]) {
    return(c(bestAlignments[1,1], bestAlignments[1,5], bestAlignments[1,6]))
  } else {
    for (i in 1:nrow(bestAlignments)) {
      if (bestAlignments[i,1] %in% row.names(wangGen) & wangGen[bestAlignments[i,1],1] == bestAlignments[i,5]) {
        return(c(bestAlignments[i,1], bestAlignments[i,5], bestAlignments[i,6]))
      }
    }
  }
})


physMap <- compact(physMap)
physMap <- data.frame(matrix(unlist(physMap), ncol = 3, byrow = T), stringsAsFactors = F)
names(physMap) <- c("Name", "Contig", "Position")
row.names(physMap) <- physMap[,1]
physMap <- physMap[,2:3]
physMap <- physMap[order(physMap[,1], as.integer(physMap[,2])),]
physMap2 <- physMap[c(which(duplicated(physMap)), which(duplicated(physMap))-1),]

## import the snp data and format
data <- read.csv("Data\\Raw\\Genotypes\\Jan_6_wheat_genotypes_curtis.csv",
                 header = T, comment.char = "", quote="", stringsAsFactors = F, row.names = 2)
colnames(data)[1:6] <- data[2,1:6]
data <- data[-2,]

data <- data[match(row.names(physMap2), row.names(data)),]

## transform the representation of the snp calls in the data
genotypes <- data[,-1:-5]
genotypes <- genotypes[,order(names(genotypes))]
genotypes <- replace(genotypes, genotypes == "c1", 0) 
genotypes <- replace(genotypes, genotypes == "C1", 0)
genotypes <- replace(genotypes, genotypes == "C2", 2)
genotypes <- replace(genotypes, genotypes == "NC", NA)

cors <- vector(length = 477)
for (i in 1:477) {
  cors[i] <- cor(as.integer(genotypes[i,]), as.integer(genotypes[i+477,]), use = "complete.obs")
}

length(which(cors > 0.99))


physMap <- physMap[-c(which(duplicated(physMap)), which(duplicated(physMap))[-which(cors > 0.99)]-1),]

save(physMap, file = "Data\\Intermediate\\Maps\\Phys\\refseqv1_best_positions.RData")
