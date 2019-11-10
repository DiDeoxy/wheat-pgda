library(SNPRelate)
library(plyr)

setwd("C:\\Users\\Max_H.DESKTOP-AJ57KB6\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## reading in the data
data <- read.csv("Data\\Raw\\Genotypes\\Jan_6_wheat_genotypes_curtis.csv",
                 header = T, comment.char = "", quote="", stringsAsFactors = F)
colnames(data)[1:6] <- data[2,1:6]

## separating the marker data from the genotype data for processing
genotypes <- data[-1:-2,-1:-6]
marker_data <- data[-1:-2,1:6]

## ordering the samples alphabetically
genotypes <- genotypes[,order(names(genotypes))]

## writing out "alphabetically" ordered cultivar names for metadata construction
# write.csv(as.vector(names(genotypes)), "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Cultivar Info\\Parsed\\metadata.csv", row.names = F)

## replaces the cluster symbols with values appropriate for the GDS format
genotypes <- replace(genotypes, genotypes == "c1", 0)
genotypes <- replace(genotypes, genotypes == "C1", 0)
genotypes <- replace(genotypes, genotypes == "C2", 2)
genotypes <- replace(genotypes, genotypes == "NC", 3)

## replaces chromosmes identifiers with numbers
marker_data$Chrom <- revalue(marker_data$Chrom, c("1A" = 1,  "1B" = 2,  "1D" = 3,  "2A" = 4,  "2B" = 5, 
                                                  "2D" = 6, "3A" = 7,  "3B" = 8,  "3D" = 9,  "4A" = 10,
                                                  "4B" = 11, "4D" = 12, "5A" = 13, "5B" = 14, "5D" = 15,
                                                  "6A" = 16, "6B" = 17, "6D" = 18, "7A" = 19, "7B" = 20,
                                                  "7D" = 21, "#N/A" = 0))

## replaces #N/A hex postion with 0 which is NA in GDSfmt
marker_data$Position <- revalue(marker_data$Position, c("#N/A" = 0))
marker_data$Position <- as.integer(as.numeric(marker_data$Position) * 1000)

## reading in metadata
metadata <- read.csv("Data\\Raw\\Parsed\\metadatav2.csv",
                     header = T, stringsAsFactors = F)

## Creates a gds format file of all wheats
samp.annot <- list(BP = metadata$Breeding.Program, Year = metadata$Date, 
                   origin = metadata$Origin, strength = metadata$Strength, colour = metadata$Colour, 
                   season = metadata$Season, designation = metadata$Designation, MC = metadata$Consensus)
snpgdsCreateGeno("Data\\Formatted\\wheat_all.gds",
                 genmat = data.matrix(genotypes),
                 sample.id = metadata$Real.Name, 
                 snp.id = as.character(marker_data$Name),
                 snp.chromosome = marker_data$Chrom,
                 snp.position = marker_data$Position,
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)