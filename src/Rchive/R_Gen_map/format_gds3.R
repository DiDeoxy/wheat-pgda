library(SNPRelate)

## reading in the data
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\snp_indices.Rdata")
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\sample_indices.Rdata")
wheat.all <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_all.gds")
sample.id <- read.gdsn(index.gdsn(wheat.all, "sample.id"))[-sample.indices]
snp.id <- read.gdsn(index.gdsn(wheat.all, "snp.id"))[snp.indices]
snp.pos <- read.gdsn(index.gdsn(wheat.all, "snp.position"))[snp.indices]
snp.chr <- read.gdsn(index.gdsn(wheat.all, "snp.chromosome"))[snp.indices]
genotypes <- read.gdsn(index.gdsn(wheat.all, "genotype"))[snp.indices, -sample.indices]
genotypes <- replace(genotypes, genotypes == 3, NA)
genotypes <- replace(genotypes, genotypes == 0, 1)
genotypes <- knncatimpute(genotypes)
genotypes <- replace(genotypes, genotypes == "NA", 3)
genotypes <- replace(genotypes, genotypes == 1, 0)
samp.annot <- read.gdsn(index.gdsn(wheat.all, "samp.annot"))
snpgdsClose(wheat.all)
for(name in names(samp.annot)) {
  samp.annot[[name]] <- factor(samp.annot[[name]][-sample.indices])
}

#writing the data out
snpgdsCreateGeno("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset_imputed.gds",
                 genmat = genotypes,
                 sample.id = sample.id, 
                 snp.id = snp.id,
                 snp.chromosome = snp.chr,
                 snp.position = snp.pos,
                 other.vars = list(samp.annot = samp.annot),
                 snpfirstdim = T)