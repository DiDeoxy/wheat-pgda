wheat <- snpgdsOpen(gdsSubset)
sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))
snp.id <- as.character(read.gdsn(index.gdsn(wheat, "snp.id")))
snp.pos <- as.integer(read.gdsn(index.gdsn(wheat, "snp.position")))
snp.chrom <- as.integer(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))

samp.annot <- read.gdsn(index.gdsn(wheat, "samp.annot"))
bp <- read.gdsn(index.gdsn(wheat, "samp.annot/BP"))
year <- read.gdsn(index.gdsn(wheat, "samp.annot/Year"))
desig <- read.gdsn(index.gdsn(wheat, "samp.annot/designation"))
mc <- read.gdsn(index.gdsn(wheat, "samp.annot/MC"))

era <- cut(as.numeric(as.character(year)),
           breaks = c(1800, 1920, 1940, 1960, 1980, 2000, 2020),
           labels = c("Pre-1920", "1921-1940", "1941-1960", "1961-1980", "1981-2000", "2001-2016")
           )
era <- addNA(era)
levels(era)[7] <- "UNKNOWN"

snpgdsClose(wheat)



