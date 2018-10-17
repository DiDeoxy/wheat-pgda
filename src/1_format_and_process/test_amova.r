genind <- read_rds(str_c(
  "Data\\Intermediate\\Adegenet\\Lr21_genind.rds"
))

marker_index <- which(names(genind$all.names) == "RAC875_c10925_1887")
table(cbind(genind$tab[, marker_index], genind$strata))

install.packages("DEMEtics")
library(DEMEtics)
data(Example.transformed)
Example.transformed