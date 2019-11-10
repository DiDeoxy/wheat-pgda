library(poppr)
library(SNPRelate)
library(ape)
library(adegenet)

base <- "C:\\Users\\Max_H.DESKTOP-AJ57KB6\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\"

## PSR 
# Per Locus data re-formatiing
load(paste(base, "Data\\Formatted\\genind_psr_hrs.RData", sep = ""))
wheat.poppr.seploc <- seploc(wheat.poppr.psr_hrs)
for (i in 1:length(wheat.poppr.seploc)){
  strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.psr_hrs)
}
# length(wheat.poppr.seploc)

wheat.amova <- list()
for (i in 1:length(wheat.poppr.seploc)) {
  if (isPoly(wheat.poppr.seploc[[i]])) {
    wheat.amova[[i]] <- poppr.amova(wheat.poppr.seploc[[i]], hier = ~psr_hrs, missing = "genotype", within = F, clonecorrect = F)
  } else {
    wheat.amova[[i]] <- NA
  }
}

# saving the amova object
save(wheat.amova, file = paste(base, "Data\\Formatted\\amova_psr_hrs.RData", sep = ""))


## HRS
# Per Locus data re-formatiing
load(paste(base, "Data\\Formatted\\genind_hrs2_hrs3.RData", sep = ""))
wheat.poppr.seploc <- seploc(wheat.poppr.hrs2_hrs3)
for (i in 1:length(wheat.poppr.seploc)){
  strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.hrs2_hrs3)
}

wheat.amova <- list()
for (i in 1:length(wheat.poppr.seploc)) {
  if (isPoly(wheat.poppr.seploc[[i]])) {
    wheat.amova[[i]] <- poppr.amova(wheat.poppr.seploc[[i]], hier = ~hrs2_hrs3, missing = "genotype", within = F, clonecorrect = F)
  } else {
    wheat.amova[[i]] <- NA
  }
}

# saving the amova object
save(wheat.amova, file = paste(base, "Data\\Formatted\\amova_hrs2_hrs3.RData", sep = ""))