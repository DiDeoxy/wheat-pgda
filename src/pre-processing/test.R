library(tidyverse)
library(SNPRelate)

# load the phys map with the genotypes
phys_map_genotypes <- read_rds(
    "Data\\Intermediate\\Maps\\phys_map_genotypes.rds")

## import categorical information on wheat varieites (market class, breeding program, year of release, phenotype, etc.)
metadata <- read_csv("Data\\Raw\\Parsed\\metadata_final.csv")

## construct the sample annotation information from the metadata
samp.annot <- list(BP = metadata$`Breeding Program`, Year = metadata$Date, 
                   origin = metadata$Origin, texture = metadata$Strength, 
                   colour = metadata$Colour, habit = metadata$Season, 
                   designation = metadata$Designation, MC = metadata$Consensus)

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