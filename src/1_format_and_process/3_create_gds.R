library(tidyverse)
library(SNPRelate)

# set the base directory
base <- file.path("data", "R")

# load the phys map with the genotypes
maps_genotypes <- read_rds(
  file.path(base, "marker_maps", "maps_genotypes.rds")
)
# select the map info
maps <- maps_genotypes %>% select(marker, chrom, phys_pos, gen_pos)
# select the genotype info and order the columns alphabetically
genotypes <- maps_genotypes %>%
  ungroup() %>%
  select(-c(marker, chrom, phys_pos, gen_pos)) %>%
  select(order(current_vars()))

# import categorical information on wheat varieites (market class,
# breeding program, year of release, phenotype, etc.)
metadata <- read_csv(
  file.path("data", "variety_info", "all_variety_info.csv") 
) %>% arrange(`Real Name`)

# # print out ordered sample names for perl cliques
# write(metadata$`Real Name`,
#   file = file.path("data", "variety_info", "ordered_names.csv"),
#   sep = "\n"
# )

# transform year into binned eras
era <- as.numeric(as.character(metadata$Date))
era <- cut(era,
  breaks = c(1800, 1920, 1940, 1960, 1980, 2000, 2020),
  labels = c(
    "Pre-1920", "1921-1940", "1941-1960", "1961-1980", "1981-2000",
    "2001-2016"
  )
)
era <- addNA(era)
levels(era)[7] <- "UNKNOWN"

## construct the sample annotation information from the metadata
samp_annot <- list(
  bp = factor(metadata$`Breeding Program`), era = factor(era),
  origin = factor(metadata$Origin), texture = factor(metadata$Texture),
  colour = factor(metadata$Colour), habit = factor(metadata$Habit),
  pheno = factor(metadata$Designation), mc = factor(metadata$Consensus)
)

# make the GDS directory if it doesnt yet exist
ifelse(
  ! dir.exists(file.path(base, "GDS")), dir.create(file.path(base, "GDS")),
  FALSE
)

## construct the SNPRelate GDS object fromt the input data with physical map
snpgdsCreateGeno(
  file.path(base, "GDS", "full_phys.gds"),
  genmat = data.matrix(genotypes),
  sample.id = metadata$`Real Name`,
  snp.id = maps$marker,
  snp.chromosome = as.integer(as.factor(maps$chrom)),
  snp.position = maps$phys_pos,
  other.vars = list(samp_annot = samp_annot),
  snpfirstdim = T
)

## construct the SNPRelate GDS object form the input data with genetic map
gen_order <- order(maps$chrom, maps$gen_pos)

snpgdsCreateGeno(
  file.path(base, "GDS", "full_gen.gds"),
  genmat = data.matrix(genotypes[gen_order, ]),
  sample.id = metadata$`Real Name`,
  snp.id = maps$marker[gen_order],
  snp.chromosome = as.integer(as.factor(maps$chrom[gen_order])),
  snp.position = maps$gen_pos[gen_order] * 100,
  other.vars = list(samp_annot = samp_annot),
  snpfirstdim = T
)
