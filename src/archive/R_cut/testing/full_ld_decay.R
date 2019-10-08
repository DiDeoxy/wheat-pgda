library(SNPRelate)
library(extrafont)
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
}
library(RColorBrewer)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")
snp.chr <- factor(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))

labels <- as.vector(t(outer(c("A", "B", "D"), c("part1", "part2"), paste, sep="_")))
labels <- as.vector(outer(as.character(1:7), labels, paste, sep=""))


count <- 1
ldMat <- list()
for(i in levels(snp.chr)) {
  ldMat[[labels[count]]] <- c(snpgdsLDMat(wheat, method = "r", snp.id = snp.id[which(snp.chr == i)], slide = 0),
                               list(pos = snp.pos[which(snp.chr == i)]))
  count = count + 1
}
snpgdsClose(wheat)

# distance <- snp.pos[which(snp.chr == 1)]
# LD.data <- ldMat[[1]]$LD
# n <- ncol(ldMat[[1]]$LD)
# HW.st<-c(C=0.5)
# HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
# tt<-summary(HW.nonlinear)
# new.rho<-tt$parameters[1]
# fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))

?nls
