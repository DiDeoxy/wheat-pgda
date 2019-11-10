# install.packages("haplotypes")
# library(haplotypes)
# install.packages("haplo.stats")
library(haplo.stats)
library(SNPRelate)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
snp.id <- as.character(read.gdsn(index.gdsn(wheat, "snp.id")))
sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))
desig <- as.factor(read.gdsn(index.gdsn(wheat, "samp.annot/designation")))
snp.chr <- as.factor(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
snpgdsClose(wheat)

index.HRS <- which(desig == "HRS")
index.SWS <- which(desig == "SWS")
index.HRS_SWS <- c(index.HRS, index.SWS)

genotypes2.3A_part2 <- geno1to2(t(genotypes[snp.chr=="14",index.HRS_SWS]))
geno_sum <- summaryGeno(genotypes2.3A_part2)
print(geno_sum)

haplo.group(desig[index.HRS_SWS][c(-9, -35, -48, -74)]=="SWS", genotypes2.3A_part2, locus.label=snp.id[snp.chr=="14"], missing.val=NA)
?haplo.group
hla.demo$label

# row.names(genotypes) <- snp.id
# colnames(genotypes) <- sample.id
# 
# genotypes <- replace(genotypes, genotypes == 0, 1)
# genotypes <- replace(genotypes, genotypes == 3, "?")
# 
# index.HRS <- which(desig == "HRS")
# index.SWS <- which(desig == "SWS")
# 
# dna.HRS <- as.dna(t(genotypes[,index.HRS]))
# dna.SWS <- as.dna(t(genotypes[,index.SWS]))
# 
# haps.HRS <- haplotype(dna.HRS, indels="missing")
# haps.SWS <- haplotype(dna.SWS, indels="missing")
# 
# pars.HRS <- parsimnet(dna.HRS, indels = "missing", prob = 0.95)
# plot(pars.HRS)
# pars.SWS <- parsimnet(dna.SWS, indels = "missing", prob = 0.80)
# plot(pars.SWS)
# 
# str(haps.SWS)
# 
# ?haplotype
# 
# ?haplotype

# nucs <- c("A", "C", "G", "T")
# 
# set.seed(1000)
# ref <- nucs[sample(4, 13192, replace = T)]
# 
# snps <- vector()
# count <- 0
# for (base in ref) {
#   if (base == "A") {
#     snps <- append(snps, sample(nucs[c(2:4)], 1))
#   } else if (base == "C") {
#     snps <- append(snps, sample(nucs[c(1,3,4)], 1))
#   } else if (base == "G") {
#     snps <- append(snps, sample(nucs[c(1,2,4)], 1))
#   } else {
#     snps <- append(snps, sample(nucs[c(1:3)], 1))
#   }
# }
# 
# both <- cbind(ref, snps)
# 
# dim(genotypes)
# 
# seqs <- vector()
# for (i in 1:dim(genotypes)[2]) {
#   seq <- vector()
#   for (j in 1:dim(genotypes)[1]) {
#     if (genotypes[j,i] == 0) {
#       seq <- append(seq, both[j,1])
#     } else if (genotypes[j,i] == 2) {
#       seq <- append(seq, both[j,2])
#     } else {
#       seq <- append(seq, NA)
#     }
#   }
#   seqs <- cbind(seqs, seq)
# }
# 
# 
# colnames(seqs) <- sample.id
# row.names(seqs) <- snp.id
# # seqs <- apply(genotypes, 2, function (seq) {
# #   lapply(seq, function (base) {
# #     if (base == 0) {
# #       seq <- append(seq, both[j,1])
# #     } else if (base == 2) {
# #       seq <- append(seq, both[j,2])
# #     } else {
# #       seq <- append(seq, NA)
# #     }
# #   }
# # })
# 
# dna <- as.dna(t(seqs))
