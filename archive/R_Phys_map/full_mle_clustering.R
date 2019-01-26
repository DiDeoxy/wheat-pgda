install.packages("mclust")
library(mclust)

setwd("C:\\Users\\Max_H.DESKTOP-AJ57KB6\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the population data
wheat <- snpgdsOpen("Data\\Formatted\\wheat_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))
snpgdsClose(wheat)

dim(genotypes)
?Mclust
mod <- Mclust(t(genotypes))
summary(mod)

png("Analysis\\Figures\\dend\\full_test.png", 
    family="Times New Roman", width = 40, height = 40, pointsize = 20, units = "in", res = 500)
plot(mod, what = "classification")
dev.off()
