library(SNPRelate)
library(plyr)
# install.packages("dbscan")
library(dbscan)
library(scrime)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the gds of the data and pullling some attributes out
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

genotypes <- replace(genotypes, genotypes == 0, 1)
genotypes <- replace(genotypes, genotypes == 3, NA)

genotypesImputed <- knncatimpute(t(genotypes))

object <- optics(genotypesImputed, minPts = 10)

for (i in 1:30) {
  print(extractDBSCAN(object, eps_cl = i))
}
dbscan <- extractDBSCAN(object, eps_cl = 28)

save(dbscan, file = "Data\\Intermediate\\dbscan\\dbscan.RData")

# library(HDclassif)
# 
# ?hddc
# model <- hddc(genotypesImputed, model = "all")
# str(model)
# model$class
