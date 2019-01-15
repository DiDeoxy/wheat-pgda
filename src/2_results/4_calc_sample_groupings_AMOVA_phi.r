library(tidyverse)
library(poppr)
library(SNPRelate)
library(ape)
library(stringr)

# genind object for amova
strata_genind <- read_rds("Data/Intermediate/Adegenet/strata_genind.rds")

# create an IBS distance object of the sample genotypes
wheat_gds <- snpgdsOpen(
  "Data/Intermediate/GDS/ld_pruned_phys_sample_subset.gds"
)
sample_dist <- as.dist(1 - snpgdsIBS(wheat_gds, autosome.only = F)$ibs)
snpgdsClose(wheat_gds)

strata <- c(
  "bp", "era", "pheno", "pheno/bp", "pheno/era", "bp/era", "era/bp", "bp_era",
  "pheno/bp_era", "clusters", "clusters/bp_era"
)

results <- tibble(
  Comparison = character(), `% Variation` = double(), `p-value` = double()
)
for (stratum in strata) {
  print(stratum)
  amova_result <- poppr.amova(
    strata_genind, hier = str_c("~", stratum) %>% as.formula(), cutoff = 0.1,
    missing = "genotype", within = FALSE, clonecorrect = FALSE, 
    dist = sample_dist, quiet = T
  )
  amova_randtest <- NULL
  if (nrow(amova_result$statphi) == 1) {
    amova_randtest <- randtest(amova_result, nrepet = 999, alter = "greater")
    results <- results %>%
      add_row(
        Comparison = stratum,
        `% Variation` = round(amova_result$componentsofcovariance$`%`[1], 2),
        `p-value` = amova_randtest$pvalue
      )
  } else {
    amova_randtest <- randtest(amova_result, nrepet = 999, output = "full") %>%
    with(
      as.krandtest(
        sim, obs, alter = c("greater", "greater"), call = call,
        names = names
      )
    )
    results <- results %>%
      add_row(
        Comparison = str_c(stratum, " (", strsplit(stratum, "/")[[1]][1], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[1], 2),
        `p-value` = amova_randtest$pvalue[3]
      )
    results <- results %>%
      add_row(
        Comparison = str_c(stratum, " (", strsplit(stratum, "/")[[1]][2], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[2], 2),
        `p-value` = amova_randtest$pvalue[2]
      )
  }
}

write_csv(results, "Results//amova_results.csv", col_names = TRUE)