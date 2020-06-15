suppressMessages(library(SNPRelate))
suppressMessages(library(randomForest))
library(plyr)

wheat.all <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_all.gds")

genotypes <- t(read.gdsn(index.gdsn(wheat.all, "genotype")))
desig <- read.gdsn(index.gdsn(wheat.all, "samp.annot/designation"))
desig <- revalue(factor(desig), c("UNKNOWN" = "N/A"))

genotypes <- replace(genotypes, genotypes == 3, NA)
genotypes <- replace(genotypes, genotypes == 0, 1)
genotypes <- as.data.frame(apply(genotypes, 2, as.factor))
genotypes <- cbind(desig, genotypes)

genotypes.imputed <- rfImpute(desig ~ ., genotypes)

genotypes.imputed <- apply(genotypes.imputed[,-1], 2, as.character)
genotypes.imputed <- apply(genotypes.imputed, 2, as.integer)
genotypes.imputed <- replace(genotypes.imputed, genotypes.imputed == 1, 0)

sample.id <- read.gdsn(index.gdsn(wheat.all, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat.all, "snp.id"))
snp.pos <- read.gdsn(index.gdsn(wheat.all, "snp.position"))
snp.chr <- read.gdsn(index.gdsn(wheat.all, "snp.chromosome"))
samp.annot <- read.gdsn(index.gdsn(wheat.all, "samp.annot"))
for(name in names(samp.annot)) {
  samp.annot[[name]] <- factor(samp.annot[[name]])
}

snpgdsClose(wheat.all)

#writing the data out
snpgdsCreateGeno("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_imputed.gds",
                 genmat = genotypes.imputed,
                 sample.id = sample.id, 
                 snp.id = snp.id,
                 snp.chromosome = snp.chr,
                 snp.position = snp.pos,
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = F)