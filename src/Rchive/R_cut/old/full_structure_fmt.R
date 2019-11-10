library(SNPRelate)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))
sample.id <- gsub(" ", "_", sample.id)
desig <- as.character(read.gdsn(index.gdsn(wheat, "samp.annot/designation")))
desig <- replace(desig, desig == "N/A", "UNKNOWN")
snpgdsClose(wheat)

#get indiv indexes
index.hrs <- which(desig == "HRS")
index.sws <- which(desig == "SWS")
index.hrw <- which(desig == "HRW")

## HRS vs SWS
index.hrs_sws <- c(index.hrs, index.sws)
genotypes.hrs_sws <- t(genotypes[,index.hrs_sws])
genotypes.hrs_sws <- cbind(as.integer(as.factor(desig[index.hrs_sws])), genotypes.hrs_sws)
colnames(genotypes.hrs_sws) <- c("", snp.id)
row.names(genotypes.hrs_sws) <- sample.id[index.hrs_sws]
write.table(genotypes.hrs_sws, "Data\\Intermediate\\Structure\\structure_hrs_sws.txt", sep = "\t", quote = F)

## HRS vs HRW
index.hrs_hrw <- c(index.hrs, index.hrw)
genotypes.hrs_hrw <- t(genotypes[,index.hrs_hrw])
genotypes.hrs_hrw <- cbind(as.integer(as.factor(desig[index.hrs_hrw])), genotypes.hrs_hrw)
colnames(genotypes.hrs_hrw) <- c("", snp.id)
rownames(genotypes.hrs_hrw) <- sample.id[index.hrs_hrw]
write.table(genotypes.hrs_hrw, "Data\\Intermediate\\Structure\\structure_hrs_hrw.txt", sep = " ", quote = F)


# ## all
# wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
# desig <- as.character(read.gdsn(index.gdsn(wheat, "samp.annot/designation")))
# desig <- replace(desig, desig == "N/A", "UNKNOWN")
# genotypes <- t(read.gdsn(index.gdsn(wheat, "genotype")))
# genotypes <- cbind(as.integer(as.factor(desig)), genotypes)
# colnames(genotypes) <- c("", snp.id)
# rownames(genotypes) <- sample.id
# snp.id <- as.character(read.gdsn(index.gdsn(wheat, "snp.id")))
# sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))
# sample.id <- gsub(" ", ".", sample.id, fixed = TRUE)
# snpgdsClose(wheat)
# write.table(genotypes, "Data\\Intermediate\\Structure\\structure.txt", sep = " ", quote = F)