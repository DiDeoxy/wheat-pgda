library(tidyverse)
library(SNPRelate)

# load the phys map with the genotypes
maps_genotypes <- read_rds(
    "Data\\Intermediate\\Maps\\maps_genotypes.rds"
)
maps <- maps_genotypes %>% select(marker, chrom, phys_pos, gen_pos)
genotypes <- maps_genotypes %>%
    ungroup() %>%
    select(-c(marker, chrom, phys_pos, gen_pos)) %>%
    select(order(current_vars()))

# import categorical information on wheat varieites (market class,
# breeding program, year of release, phenotype, etc.)
metadata <- read_csv("Data\\Raw\\Parsed\\metadata_final.csv") %>%
    arrange(`Real Name`)

# print out ordered sample names for perl cliques
write(metadata$`Real Name`,
    file = "Data\\Intermediate\\ordered_names.txt",
    sep = "\n"
)

## construct the sample annotation information from the metadata
samp_annot <- list(
    BP = metadata$`Breeding Program`, Year = metadata$Date,
    origin = metadata$Origin, texture = metadata$Texture,
    colour = metadata$Colour, habit = metadata$Habit,
    designation = metadata$Designation, MC = metadata$Consensus
)

## construct the SNPRelate GDS object fromt the input data with physical map
snpgdsCreateGeno("Data\\Intermediate\\GDS\\full_phys.gds",
    genmat = data.matrix(genotypes),
    sample.id = metadata$`Real Name`,
    snp.id = maps$marker,
    snp.chromosome = as.integer(as.factor(maps$chrom)),
    snp.position = maps$phys_pos,
    other.vars = list(samp_annot = samp_annot),
    snpfirstdim = T
)

## construct the SNPRelate GDS object form the input data with genetic map
snpgdsCreateGeno("Data\\Intermediate\\GDS\\full_gen.gds",
    genmat = data.matrix(genotypes),
    sample.id = metadata$`Real Name`,
    snp.id = maps$marker,
    snp.chromosome = as.integer(as.factor(maps$chrom)),
    snp.position = maps$gen_pos,
    other.vars = list(samp_annot = samp_annot),
    snpfirstdim = T
)
