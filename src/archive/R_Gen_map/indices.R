library(SNPRelate)

wheat.all <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_all.gds")
snp.id <- read.gdsn(index.gdsn(wheat.all, "snp.id"))
sample.id <- read.gdsn(index.gdsn(wheat.all, "sample.id"))
genotypes <- read.gdsn(index.gdsn(wheat.all, "genotype"))

## SNP indices
set.seed(1000)
snp.set <- snpgdsLDpruning(wheat.all, ld.threshold=0.3, autosome.only = F,
                           maf = 0.05, missing.rate = 0.1, slide.max.bp = 345)

snp.set.id <- unlist(snp.set)
snp.indices <- match(snp.set.id, snp.id)
save(snp.indices, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\snp_indices.Rdata")

## Sample indices
IBS <- snpgdsIBS(wheat.all)
snpgdsClose(wheat.all)

pairs <- which(IBS$ibs >= 0.99, arr.ind = T)
indices <- vector()
for (i in 1:dim(pairs)[1]) {
  if (pairs[i,1] == pairs[i,2]){
    indices <<- c(indices, i)
  }
}
pairs <- pairs[-indices,]
pairs

## send pairs to cliques.pl, remove near duplicates
# NILs <- c("PT434","BW811","Conway","AC Minto 1","Avocet 1","BW275 1","PT616","BW427 1","BW942","BW922",
#           "BW948","Carberry 1","CDC Stanley 1","RL4452","Winalta","PT754","SWS349","Somerset 1","Stettler 1","AC Andrew",
#           "SWS363","CDC Makwa","Neepawa","AC Reed 1","SWS87","SWS345","SWS363","SWS390","SWS408", "SWS410")

NILs <- c("AC Elsa", "AC Minto 1", "AC Reed 1", "SWS87", "Avocet 1", "BW275 1", "PT616", "BW427 1", "BW942", "BW948", "Carberry 1", "Neepawa", 
          "CDC Stanley 1", "RL4452", "Winalta", "Somerset 1", "Stettler 1", "SWS363", "SWS241", "SWS345", "SWS363", "SWS390", "SWS408", "SWS410")

sample.indices <- match(NILs, sample.id)
save(sample.indices, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\sample_indices.Rdata")

## making the default hclust plot
plot(hclust(as.dist(1-cor(genotypes[snp.indices, -sample.indices], use = "pairwise.complete.obs")), method = "average"))
