library(poppr)
library(SNPRelate)
library(ape)
library(stringr)


## loading the gds of the data and pullling some attributes out
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("src\\R_functions\\data_loading.R")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
distances <- 1 - snpgdsIBS(wheat, autosome.only = F)$ibs
snpgdsClose(wheat)

genind <- get(load(paste0("Data\\Intermediate\\Adegenet\\genind_wheat.RData")))

results <- list()
for (grouping in c(
  "HS", "RW", "SW",
  "RWHS/SW", "SW/RWHS",
  "SWHS/RW", "RW/SWHS",
  "SWRW/HS", "HS/SWRW",
  "desig", "era", "bp", "HRSOther",
  "SWHS/RW/bp",
  "hdbscan2", "hdbscan6", "hdbscan2/hdbscan6",
  "dend"
)[1:3]) {
  print(grouping)
  amovaResult <- poppr.amova(genind,
    hier = as.formula(paste0("~", grouping)), cutoff = 0.1,
    missing = "genotype", within = F, clonecorrect = F,
    dist = as.dist(distances), quiet = T
  )
  print(amovaResult)
  if (grepl("/", grouping)) {
    amovaRandtest <- randtest(amovaResult, nrepet = 999, output = "full")
    amovaRandtest <- with(amovaRandtest, as.krandtest(sim, obs, alter = c("two-sided", "two-sided", "two-sided"), call = call, names = names))
    results[[grouping]][[1]] <- cbind(amovaResult$componentsofcovariance, "direction" = rev(c("NA", amovaRandtest$alter)), "p-value" = rev(c("NA", amovaRandtest$pvalue)))
  } else {
    amovaRandtest <- randtest(amovaResult, nrepet = 999, alter = "greater")
    results[[grouping]][[1]] <- cbind(amovaResult$componentsofcovariance, "direction" = rev(c("NA", "NA", amovaRandtest$alter)), "p-value" = rev(c("NA", "NA", amovaRandtest$pvalue)))
  }
}
results

lapply(results, function(x) {
  write.table(data.frame(x), "amova_results.csv", append = T, sep = ",")
})
