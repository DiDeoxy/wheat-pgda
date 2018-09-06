wheat <- snpgdsOpen(gds)
sample_id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))
snp_id <- as.character(read.gdsn(index.gdsn(wheat, "snp.id")))
snp_pos <- as.integer(read.gdsn(index.gdsn(wheat, "snp.position")))
snp_chrom <- as.integer(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))

samp_annot <- read.gdsn(index.gdsn(wheat, "samp_annot"))
bp <- read.gdsn(index.gdsn(wheat, "samp_annot/BP"))
year <- read.gdsn(index.gdsn(wheat, "samp_annot/Year"))
desig <- read.gdsn(index.gdsn(wheat, "samp_annot/designation"))
mc <- read.gdsn(index.gdsn(wheat, "samp_annot/MC"))

era <- cut(as.numeric(as.character(year)),
           breaks = c(1800, 1920, 1940, 1960, 1980, 2000, 2020),
           labels = c("Pre-1920", "1921-1940", "1941-1960", "1961-1980", "1981-2000", "2001-2016")
           )
era <- addNA(era)
levels(era)[7] <- "UNKNOWN"

snpgdsClose(wheat)



