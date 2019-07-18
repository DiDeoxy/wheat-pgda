source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(pgda, "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_point", "labs",  "scale_colour_gradientn", 
  "scale_size_continuous", "theme", "unit", "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(readr, "read_rds")
import::from(Rfast, "rowMaxs")
import::from(scrime, "rowTables")
import::from(stringr, "str_c")
import::from(tibble, "as_tibble", "tibble")

# load the data from the gds object
phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

# recode the genotypes and count
genos <- replace(phys_data$genotypes, phys_data$geno == 3, NA)

# make a tibble with the relevant data
coding <- c(0, 2)
allele_counts <- rowTables(genos, coding)
snp_data <- tibble(
  chrom = phys_data$snp$chrom,
  phys_pos_mb = phys_data$snp$pos / 1e6,
  gen_pos_cm = gen_data$snp$pos[match(phys_data$snp$id, gen_data$snp$id)] / 100,
  josts_d := read_rds(josts_d),
  overall_mjaf =
    rowMaxs(allele_counts, value = TRUE) / rowSums(allele_counts),
)

# get the cluster groups
cluster <- read_rds(hdbscan)$cluster
subs_genos <- list(
  chrs_mjaf = genos[,
    which(phys_data$sample$annot$pheno == "HRS" & cluster == 5)
  ],
  chrw_mjaf = genos[,
    which(phys_data$sample$annot$pheno == "HRW" & cluster == 1)
  ],
  csws_mjaf = genos[,
    which(phys_data$sample$annot$pheno == "SWS" & cluster == 2)
  ]
)

mja <- rowMaxs(rowTables(genos, coding))
mjafs_by_pop <- lapply(subs_genos, function (sub_genos) {
  genos_counts <- rowTables(sub_genos, coding)
  max_genos <- genos_counts[cbind(seq_along(mja), mja)]
  max_genos / rowSums(genos_counts)
}) %>% do.call(cbind, .)

# add the genes positons to the regions table
snp_data <- snp_data %>%
  cbind(mjafs_by_pop %>% round(4)) %>%
  split(.$chrom)

# create a function for making a gradient of colours
colour_gradient <- colorRampPalette(colour_set[c(1, 5, 3, 2, 4)])(100)

# calc the lengths of the different genomes and homoeologous sets
max_phys_lengths <- span_by_chrom(
  phys_data$snp$chrom, phys_data$snp$pos
) %>% max_lengths() / 1e6
max_gen_lengths <- span_by_chrom(
  gen_data$snp$chrom, gen_data$snp$pos
) %>% max_lengths() / 100


lapply(snp_data, function (chrom) {
  
})