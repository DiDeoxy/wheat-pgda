library(tidyverse)
library(poppr)
library(SNPRelate)
library(ape)
library(stringr)


## loading the gds of the data and pullling some attributes out
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds"
source("src\\R_functions\\data_loading.R")

full <- snpgdsOpen(
  "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds"
)
full_dist <- 1 - snpgdsIBS(full, autosome.only = F)$ibs
snpgdsClose(full)

full_genind <- read_rds("Data\\Intermediate\\Adegenet\\full_genind.rds")

listing <- function(reps, grouping, amova_result, amova_randtest) {
  results[[grouping]][[1]] <- cbind(amova_result$componentsofcovariance,
    "direction" = rev(c(
      rep("NA", reps),
      amova_randtest$alter
    )),
    "p-value" = rev(c(
      rep("NA", reps),
      amova_randtest$pvalue
    ))
  )
  return(results)
}

results <- list()
for (grouping in c(
  "era", "bp", "texture", "colour", "habit", "phenotype", "hrs", "sws", "hrw",
  "clusters", "cluster1", "cluster2", "cluster3", "cluster4", "cluster5"
)) {
  print(grouping)
  amova_result <- poppr.amova(full_genind,
    hier = as.formula(paste0(
      "~", grouping
    )),
    cutoff = 0.1, missing = "genotype", within = F,
    clonecorrect = F, dist = as.dist(full_dist),
    quiet = T
  )
  print(amova_result)
  if (grepl("/", grouping)) {
    amova_randtest <- randtest(amova_result, nrepet = 999, output = "full")
    amova_randtest <- with(
      amova_randtest,
      as.krandtest(sim, obs,
        alter = c(
          "two-sided", "two-sided", "two-sided"
        ),
        call = call, names = names
      )
    )
    results <- listing(1, grouping, amova_result, amova_randtest)
  } else {
    amova_randtest <- randtest(amova_result, nrepet = 999, alter = "greater")
    results <- listing(2, grouping, amova_result, amova_randtest)
  }
}
results

for (name in names(results)) {
  write_csv(data.frame(results[name]),
    "Results//amova_results.csv", append = TRUE, col_names = T
  )
}