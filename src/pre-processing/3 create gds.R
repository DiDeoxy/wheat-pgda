library(SNPRelate)
library(plyr)
library(dplyr)

## importing the snp-chip genotype calls from a CSV file
data <- read.csv("Data\\Raw\\Genotypes\\Jan_6_wheat_genotypes_curtis.csv",
                 header = T, comment.char = "", quote="", stringsAsFactors = F, row.names = 2)
## formatting the data frame
colnames(data)[1:6] <- data[2,1:6]
data <- data[-2,]

## import the constructed physical map, from a previous step not shown
## used alignment of snp probes to new wheat reference genome and genetic
## map linkage info to construct a physical map
load("Data\\Intermediate\\Maps\\Phys\\refseqv1_best_positions.R")

## subset and order the data rows by the phys map order
data <- data[match(row.names(physMap), row.names(data)),]

## transform the representation of the snp calls in the data
genotypes <- data[,-1:-5]
genotypes <- genotypes[,order(names(genotypes))]
genotypes <- replace(genotypes, genotypes == "c1", 0) 
genotypes <- replace(genotypes, genotypes == "C1", 0)
genotypes <- replace(genotypes, genotypes == "C2", 2)
genotypes <- replace(genotypes, genotypes == "NC", NA)

## import categorical information on wheat varieites (market class, breeding program, year of release, phenotype, etc.)
metadata <- read.csv("Data\\Raw\\Parsed\\metadatav2.csv",
                     header = T, stringsAsFactors = T)

## construct the sample annotation information from the metadata
samp.annot <- list(BP = metadata$Breeding.Program, Year = metadata$Date, 
                   origin = metadata$Origin, strength = metadata$Strength, colour = metadata$Colour, 
                   season = metadata$Season, designation = metadata$Designation, MC = metadata$Consensus)

## construct the SNPRelate GDS object fromt the input data with physical map
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_phys.gds",
                 genmat = data.matrix(genotypes),
                 sample.id = metadata$Real.Name, 
                 snp.id = row.names(physMap),
                 snp.chromosome = as.integer(as.factor(physMap$Contig)),
                 snp.position = physMap$Position,
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)

## construct the SNPRelate GDS object form the input data with genetic map
snpgdsCreateGeno("Data\\Intermediate\\GDS\\wheat_gen.gds",
                 genmat = data.matrix(genotypes),
                 sample.id = metadata$Real.Name, 
                 snp.id = row.names(physMap),
                 snp.chromosome = as.integer(as.factor(physMap$Contig)),
                 snp.position = data$Position,
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)