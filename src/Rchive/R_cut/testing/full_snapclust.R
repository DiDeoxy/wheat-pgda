library(adegenet)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

load("Data\\Intermediate\\Adegenet\\genind_all.RData")
load("Data\\Intermediate\\CCP\\full_kmeans.RDATA")

sc.wheat <- snapclust(wheat.poppr.all, 6, pop.ini = "kmeans")

save(sc.wheat, file = "Data\\Intermediate\\full_snapclust.Rdata")
