library(SNPRelate)
library(plyr)
library(dendextend)
# install.packages("erer")
library(erer)
library(fpc)

pamk <- pamk(dist(1-abs(cor(genotypes[kept.indices,]))), krange=2:15, diss = T)

## loading the gds of the data and pullling some attributes out
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data.gds")

# Removing SNPs with high MR and low MAF
snpset <- snpgdsSelectSNP(wheat.data, autosome.only = F, maf = 0.05, missing.rate = 0.1)
snpset.id <- unlist(snpset)
informative <- match(snpset.id, snp.id)
genotypes <- read.gdsn(index.gdsn(wheat.data, "genotype"))[informative,]

counts <- list()
counts[[1]] <- "[loci]=24841"
counts[[2]] <- "[populations]=2"
count <- 3
by(t(genotypes), factor(pamk$pamobject$clustering), function (x) {
  counts[[count]] <<- data.frame()
  for(i in 1:length(x)) {
    counts[[count]][i,1] <<- length(x[,i])
    counts[[count]][i,2] <<- sum(x[,i] == 0)
    counts[[count]][i,3] <<- sum(x[,i] == 2)
    counts[[count]][i,4] <<- sum(x[,i] == 3)
  }
  count <<- count + 1
})

sink("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\pop_snps.txt")
options(max.print=1000000)
counts
sink()