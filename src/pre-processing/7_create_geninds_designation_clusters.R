library(SNPRelate)
library(adegenet)
library(dbscan)

# ## loading the gds of the data and pullling some attributes out
# gds <- "Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds"
# source("src\\functions\\data_loading.R")
# 
# load("Data\\Intermediate\\dbscan\\wheatHdbscan.RData")
# wheatHdbscan6 <- factor(wheatHdbscan$cluster)
# wheatHdbscan2 <- as.integer(wheatHdbscan6)
# wheatHdbscan2[which(wheatHdbscan2 != 6)] <- 1
# 
# ## for the major three groups
# indexNotCanadian <- c(which(dbscan6 != 1))
# indexHrsSws <- which(desig == "HRS")
# indexHrsSws <- indexHrsSws[-which(indexHrsSws %in% indexNotCanadian)]
# indexHrsSws <- c(indexHrsSws, which(desig == "SWS"))
# indexHrsHrw <- which(desig == "HRS")
# indexHrsHrw <- indexHrsHrw[-which(indexHrsHrw %in% indexNotCanadian)]
# indexHrsHrw <- c(indexHrsHrw, which(desig == "HRW"))
# indexHrwSws <- c(which(desig == "HRW"), which(desig == "SWS"))
# 
# grouping <- list(list(desig, desig, desig),
#                  list(indexHrsSws, indexHrsHrw, indexHrwSws),
#                  list("HRS_SWS", "HRS_HRW", "HRW_SWS"))
# 
# for (i in 1:length(grouping)) {
#   type <- grouping[[1]][[i]]
#   index <- grouping[[2]][[i]]
#   name <- grouping[[3]][[i]]
#   
#   strata <- data.frame(type[index])
#   colnames(strata) <- name
#   
#   genind <- df2genind(t(data.frame(genotypes[,index])),
#                       ind.names = as.character(sample.id)[index],
#                       loc.names = snp.id, NA.char = "N", ploidy = 1, type = "codom",
#                       ncode = 1, strata = strata)
#   
#   save(genind, file = paste0("Data\\Intermediate\\Adegenet\\genind_", name, ".RData"))
# }
# 
# ## resistance genes
# genePres <- read.csv("Data\\Intermediate\\Gene_presence_randhawa.csv", header = F, row.names = 1, stringsAsFactors = F)
# for (gene in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
#   indivs <- c(which(genePres$V2 == gene), which(genePres$V3 == gene), which(genePres$V4 == gene))
#   indexIndivs <- match(rownames(genePres[indivs, ]), sample.id)
#   indexIndivs <- indexIndivs[!is.na(indexIndivs)]
#   indexNotIndivs <- match(rownames(genePres[-indivs, ]), sample.id)
#   indexNotIndivs <- indexNotIndivs[!is.na(indexNotIndivs)]
# 
#   strata <- data.frame(c(rep(gene, length(indexIndivs)), rep(paste0("not", gene), length(indexNotIndivs))))
#   colnames(strata) <- gene
# 
#   genind <- df2genind(t(data.frame(genotypes[,c(indexIndivs, indexNotIndivs)])),
#                       ind.names = as.character(sample.id)[c(indexIndivs, indexNotIndivs)],
#                       loc.names = snp.id, NA.char = "N", ploidy = 1, type = "codom",
#                       ncode = 1, strata = strata)
# 
#   save(genind, file = paste0("Data\\Intermediate\\Adegenet\\genind_", gene, ".RData"))
# }