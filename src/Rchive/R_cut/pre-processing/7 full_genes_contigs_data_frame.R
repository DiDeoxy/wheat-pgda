library(SNPRelate)
library(plyr)

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds"
source("src\\functions\\data_loading.R")

makeGenesContigsDataframe <- function (fileIn, fileOut) {
  genes <- read.csv(fileIn, header = F, stringsAsFactors = F, col.names = c("Name", "ID", "Contig", "Pos", "sleng", "aleng", "%id"))[,c(3,4,1)]
  genes <- genes[order(genes[,1]),]
  row.names(genes) <- NULL
  
  chromes <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
  part2 <- paste0(chromes, "_part2")
  part2Start <- cbind(part2, c(471304005, 438720154, 452179604, 462376173, 453218924, 462216879, 454103970,
                               448155269, 476235359, 452555092, 451014251, 451004620, 453230519, 451372872,
                               451901030, 452440856, 452077197, 450509124, 450046986, 453822637, 453812268))
  for (i in 1:nrow(part2Start)) {
    genes[which(genes$Contig == part2Start[i, 1]), 2] = genes[which(genes$Contig == part2Start[i, 1]), 2] + as.integer(part2Start[i, 2])
  }
  for (row in 1:nrow(genes)) {
    genes$Contig[row] = substr(genes$Contig[row], 1, 2)
  }
  
  replace <- 1:21
  names(replace) <- chromes
  genes$Contig <- revalue(genes$Contig, replace)
  
  genesContigs <- data.frame()
  count <- 1
  number <- 0
  for (i in 1:length(snp.chrom)) {
    if (count == snp.chrom[i]) {
      hits <- genes[which(genes$Contig == count),]
      number <- nrow(hits)
      if (nrow(hits) == 0) {
        temp <- data.frame(NA, NA, NA)
        names(temp) <- names(genes)
        genesContigs <- rbind(genesContigs, temp)
      } else {
        for (j in 0:dim(hits)[1]) {
          genesContigs <- rbind(genesContigs, hits[j,])
        }
      }
      if (number > 0) { number <- number - 1 }
      count = count + 1
    } else {
      if (number == 0) {
        temp <- data.frame(NA, NA, NA)
        names(temp) <- names(genes)
        genesContigs <- rbind(genesContigs, temp)
      } else {
        number <- number - 1
      }
    }
  }
  save(genesContigs, file = fileOut)
}

makeGenesContigsDataframe("Data\\Intermediate\\Aligned_genes\\top_main_genes_selected.csv", 
                          "Data\\Intermediate\\Aligned_genes\\top_main_genes_contigs.RData")
makeGenesContigsDataframe("Data\\Intermediate\\Aligned_genes\\top_resistance_genes_selected.csv",
                          "Data\\Intermediate\\Aligned_genes\\top_resistance_genes_contigs.RData")
makeGenesContigsDataframe("Data\\Intermediate\\Aligned_genes\\known_genes_groups.csv",
                          "Data\\Intermediate\\Aligned_genes\\known_genes_groups.RData")
