library(tidyverse)
library(poppr)
library(SNPRelate)
library(ape)
library(stringr)

# genind object for amova
strata_genind <- read_rds("Data/Intermediate/Adegenet/strata_genind.rds")

# create an IBS distance object of the sample genotypes
wheat_gds <- snpgdsOpen(
  "Data/Intermediate/GDS/phys_subset_sample_ld_pruned.gds"
)
sample_dist <- as.dist(1 - snpgdsIBS(wheat_gds, autosome.only = F)$ibs)
snpgdsClose(wheat_gds)

results <- tibble()
for (stratum in names(strata_genind$strata)) {
  print(stratum)
  amova_result <- poppr.amova(
    strata_genind, hier = str_c("~", stratum) %>% as.formula(), cutoff = 0.1,
    missing = "genotype", within = FALSE, clonecorrect = FALSE, 
    dist = sample_dist, quiet = T
  )
  amova_randtest <- randtest(amova_result, nrepet = 999, alter = "greater")
  results <- c(
    stratum, round(amova_result$componentsofcovariance$`%`[1], 2),
    amova_randtest$pvalue
    ) %>%
    rbind(results, ., stringsAsFactors = FALSE)
}
names(results) <- c("Comparison", "Phi", "p-value")
results

write_csv(results, "Results//amova_results.csv", col_names = TRUE)
