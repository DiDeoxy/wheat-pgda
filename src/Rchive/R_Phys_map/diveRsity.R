library(SNPRelate)
# install.packages("diveRsity")
library(diveRsity)

wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_imputed_sample_subset.gds")
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.chr <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
samp.annot <- read.gdsn(index.gdsn(wheat, "samp.annot"))
snpgdsClose(wheat)

genotypes <- replace(genotypes, genotypes == 0, "GG")
genotypes <- replace(genotypes, genotypes == 2, "TT")
genotypes <- replace(genotypes, genotypes == 3, "--")
genotypes <- rbind(snp.id, genotypes)
genotypes <- cbind(c("snp", sample.id), genotypes)

snp2gen(genotypes, prefix_length = 0)

genotypes[1:10,1:10]
dim(genotypes)
length(snp.id)
dim(genotypes[-1,])
length(sample.id)

div <- fastDivPart("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\snp2gen_converted.gen")
head(div$standard)
