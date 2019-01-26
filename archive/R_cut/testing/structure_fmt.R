## outputs data in a format readable by STRUCTURE
library(SNPRelate)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

markers <- rbind(snp.id)
data <- cbind(sample.id, as.numeric(desig), t(data.matrix(genotypes)))

write.table(markers, "Data\\Intermediate\\Structure\\structure_both.txt",
            col.names = F, row.names = F, sep = " ", quote = F)
write.table(data, "Data\\Intermediate\\Structure\\structure_both.txt",
            col.names = F, row.names = F, sep = " ", quote = F, append = T)

