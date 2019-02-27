suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(adegenet))
suppressPackageStartupMessages(library(plyr))
source(file.path("src", "file_paths.R"))
source(parse_create_gds)

# create geninds for calclulating Josts D values
wheat_data <- parse_gds(phys_gds)

# making the genotype data palatable by genind
wheat_data$genotypes <- wheat_data$genotypes %>%
  replace(. == 0, "A") %>%
  replace(. == 2, "B") %>%
  replace(. == 3, "")

cluster <- read_rds(hdbscan)$cluster

# making indices in order to create subsets containing just the varieties of
# the phenotype/cluster groups we are interested in
index_chrs <- which(
  wheat_data$sample$annot$pheno == "HRS" & cluster == 5
)
index_chrw <- which(
  wheat_data$sample$annot$pheno == "HRW" & cluster == 1
)
index_csws <- which(
  wheat_data$sample$annot$pheno == "SWS" & cluster == 2
)

index_chrs_chrw_csws <- c(index_chrs, index_chrw, index_csws)

set.seed(100)
sampled_index_chrs <- index_chrs[length(index_chrs) %>% sample(31)]
set.seed(100)
sampled_index_csws <- index_csws[length(index_csws) %>% sample(31)]
index_sampled_chrs_chrw_csws <- c(sampled_index_chrs, index_chrw, sampled_index_csws)

grouping <- list(
  index = list(index_chrs_chrw_csws, index_sampled_chrs_chrw_csws),
  name = c("chrs_chrw_csws", "sampled_chrs_chrw_csws")
)

for (i in 1:length(grouping)) {
  index <- grouping$index[[i]]
  name <- grouping$name[i]
  print(name)
  pop <- as.character(wheat_data$sample$annot$pheno[index])

  subset_genind <- df2genind(t(data.frame(wheat_data$genotypes[, index])),
    ind.names = wheat_data$sample$id[index], loc.names = wheat_data$snp$id,
    ploidy = 1, type = "codom", ncode = 1,
    pop = pop
  )

  write_rds(subset_genind, path = file.path(geninds, str_c(name, ".rds")))
}

################################################################################
# Make a geneind with all strata to perform AMOVA analyses on them
wheat_data <- parse_gds(ld_gds)

cluster <- read_rds(hdbscan)$cluster

# making the genotype data palatable by genind
wheat_data$genotypes <- wheat_data$genotypes %>%
  replace(. == 0, "A") %>%
  replace(. == 2, "B") %>%
  replace(. == 3, "")

# make unknown and small groups uniformaly NA
wheat_data$sample$annot$bp <- revalue(
  wheat_data$sample$annot$bp,
  replace = c("FOREIGN" = "N/A", "USA" = "N/A")
)

# create vectors giving membership or not of varieties to each of the three
# largest phenotypic groups
major_phenos <- list()
for (pheno in c("HRS", "HRW", "SWS")) {
  major_phenos[[pheno]] <- as.character(wheat_data$sample$annot$pheno)
  subset_index <- which(wheat_data$sample$annot$pheno == pheno)
  major_phenos[[pheno]][-subset_index] <- "Other"
}

# create membership vectors of each cluster vs all other
cluster_isolated <- list()
for (cluster_num in c(1, 2, 3, 4, 5)) {
  cluster_isolated[[cluster_num]] <- cluster
  cluster_isolated[[cluster_num]][-which(cluster == cluster_num)] <- 0
}

# create a data frame of the different strata
strata <- tibble(
  era = wheat_data$sample$annot$era,
  bp = wheat_data$sample$annot$bp,
  bp_era = str_c(wheat_data$sample$annot$era, wheat_data$sample$annot$bp),
  texture = wheat_data$sample$annot$texture,
  colour = wheat_data$sample$annot$colour,
  habit = wheat_data$sample$annot$habit,
  pheno = wheat_data$sample$annot$pheno,
  hrs = as.factor(major_phenos$HRS),
  hrw = as.factor(major_phenos$HRW),
  sws = as.factor(major_phenos$SWS),
  clusters = as.factor(cluster),
  cluster1 = as.factor(cluster_isolated[[1]]),
  cluster2 = as.factor(cluster_isolated[[2]]),
  cluster3 = as.factor(cluster_isolated[[3]]),
  cluster4 = as.factor(cluster_isolated[[4]]),
  cluster5 = as.factor(cluster_isolated[[5]])
)

df2genind(
  t(wheat_data$genotypes),
  ind.names = wheat_data$sample$id,
  loc.names = wheat_data$snp$id, ploidy = 1, type = "codom", ncode = 1,
  strata = strata, pop = as.factor(cluster)
) %>%
  write_rds(path = str_c(file.path(geninds, "strata.rds")))