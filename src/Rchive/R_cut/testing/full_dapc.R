library(adegenet)
library(SNPRelate)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA")

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

load("Data\\Intermediate\\Adegenet\\genind_wheat.RData")

grp <- find.clusters(genind, max.n.clust = 40, n.pca = 500, n.clust = 4)


dapc1 <- dapc(genind, grp$grp, n.pca = 100)
scatter(dapc1, grp = desig)
assignplot(dapc1)
dapc1$grp
?scatter.dapc
loadingplot(dapc1$var.contr, axis = 1, thresh = 0.0015)
genind$all.names
which(names(genind$all.names) == "wsnp_Ex_c11085_17973016")
?df2genind

sizes = c(6,6)
pops = sum(sizes)
perms = factorial(pops) ???
prod(sapply(sizes,FUN = factorial)) ???
prod(sapply(table(sizes),FUN=factorial))
data.frame(perms,min.p=1/perms)