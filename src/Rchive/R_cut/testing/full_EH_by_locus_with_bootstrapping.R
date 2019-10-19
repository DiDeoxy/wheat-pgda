library(SNPRelate)
library(extrafont)
library(ape)
# install.packages("pastecs")
# library(pastecs)
# install.packages("zoo")
# library(zoo)
library(boot)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## Functions
EH <- function (data, indices) {
  apply(t(data[indices,]), 1, function (row) {
    2*((sum(row == 0)/sum(row == 0 | 2)) * (sum(row == 2)/sum(row == 0 | 2)))
  })
}

signif <- function (x) { 
  if (x[1] < x[2] | x[1] > x[3]) { return(TRUE) } else { return(FALSE) }
}

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("Data\\Formatted\\wheat_phys_subset_sample.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.chrom <- factor(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
pop.code.de <- read.gdsn(index.gdsn(wheat, "samp.annot/designation"))
snpgdsClose(wheat)

##
index.hrs <- which(pop.code.de == "HRS")
index.sws <- which(pop.code.de == "SWS")
index.hrw <- which(pop.code.de == "HRW")

##
eh.boot.hrs <- boot(data = t(genotypes[,index.hrs]), statistic = EH, R = 100)
eh.boot.hrw <- boot(data = t(genotypes[,index.hrw]), statistic = EH, R = 100)
eh.boot.sws <- boot(data = t(genotypes[,index.sws]), statistic = EH, R = 100)

##
difs.eh.boot.hrs_hrw <- eh.boot.hrs$t - eh.boot.hrw$t
difs.eh.hrs_hrw <- eh.boot.hrs$t0 - eh.boot.hrw$t0
difs.eh.boot.hrs_sws <- eh.boot.hrs$t - eh.boot.sws$t
difs.eh.hrs_sws <- eh.boot.hrs$t0 - eh.boot.sws$t0

##
quants.hrs_hrw <- apply(difs.eh.boot.hrs_hrw, 2, function (x) { quantile(x, probs = c(0.05, 0.95)) })
quants.hrs_sws <- apply(difs.eh.boot.hrs_sws, 2, function (x) { quantile(x, probs = c(0.05, 0.95)) })

##
signif.difs.hrs_hrw <- apply(rbind(difs.eh.hrs_hrw, quants.hrs_hrw), 2, signif)
which(signif.difs.hrs_hrw)

signif.difs.hrs_hrw <- apply(rbind(difs.eh.hrs_sws, quants.hrs_sws), 2,  signif)
which(signif.difs.hrs_hrw)

