suppressMessages(library(SNPRelate))
library(plyr)
suppressMessages(library(dplyr))

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## import the snp data and format
data <- read.csv("Data\\Raw\\Genotypes\\Jan_6_wheat_genotypes_curtis.csv",
                 header = T, comment.char = "", quote="", stringsAsFactors = F, row.names = 2)
colnames(data)[1:6] <- data[2,1:6]
data <- data[-2,]

## import the constructed physical map and remove duplicates
physMap <- read.csv("Data\\Intermediate\\Maps\\Phys\\pozniak_physical_map168.csv",
                     header = F, col.names = c("Name", "Contig", "Position"), row.names = 1, stringsAsFactors = F)
physMap <- physMap[order(physMap$Contig, physMap$Position),]
physMap <- physMap[-c(which(duplicated(physMap)), which(duplicated(physMap))-1),]

## combine parts into single chromosome
chromes <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
part2 <- as.vector(t(outer(chromes, "part2", paste, sep="_")))
part2Start <- cbind(part2, c(471304005, 438720154, 452179604, 462376173, 453218924, 462216879, 454103970,
                             448155269, 476235359, 452555092, 451014251, 451004620, 453230519, 451372872, 
                             451901030, 452440856, 452077197, 450509124, 450046986, 453822637, 453812268))
for (i in 1:nrow(part2Start)) {
  physMap[which(physMap$Contig == part2Start[i, 1]), 2] = physMap[which(physMap$Contig == part2Start[i, 1]), 2] + as.integer(part2Start[i, 2])
}
for (row in 1:nrow(physMap)) {
  physMap$Contig[row] = substr(physMap$Contig[row], 1, 2)
}

## subset and order the data rows by the phys map order
data <- data[match(row.names(physMap), row.names(data)),]

## transform the representation of the snp calls in the data
genotypes <- data[,-1:-5]
genotypes <- genotypes[,order(names(genotypes))]
genotypes <- replace(genotypes, genotypes == "c1", 0)
genotypes <- replace(genotypes, genotypes == "C1", 0)
genotypes <- replace(genotypes, genotypes == "C2", 2)
genotypes <- replace(genotypes, genotypes == "NC", 3)

## import metadata
metadata <- read.csv("Data\\Raw\\Parsed\\metadatav2.csv",
                     header = T, stringsAsFactors = T)

## construct the sample annotation information from the metadata
samp.annot <- list(BP = metadata$Breeding.Program, Year = metadata$Date, 
                   origin = metadata$Origin, strength = metadata$Strength, colour = metadata$Colour, 
                   season = metadata$Season, designation = metadata$Designation, MC = metadata$Consensus)

## construct the GDS object fromt the input data
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_phys.gds",
                 genmat = data.matrix(genotypes),
                 sample.id = metadata$Real.Name, 
                 snp.id = row.names(physMap),
                 snp.chromosome = as.integer(as.factor(physMap$Contig)),
                 snp.position = physMap$Position,
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)