library(SNPRelate)
library(plyr)
library(adegenet)
# library(dbscan)

## loading the gds of the data and pullling some attributes out
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("src\\functions\\data_loading.R")

# making the data palatable by genind
genotypes <- replace(genotypes, genotypes == 0, "0")
genotypes <- replace(genotypes, genotypes == 2, "2")
genotypes <- replace(genotypes, genotypes == 3, "N")

## for making subsets of the samples for the different categories excluding unknown samples
indexMC <- c(which(mc == "CERS"), which(mc == "CPSR"), which(mc == "CPSW"), which(mc == "CWES"), 
             which(mc == "CWGP"), which(mc == "CWHWS"), which(mc == "CWRS"), which(mc == "CWRW"), 
             which(mc == "CWSWS"), which(mc == "CXHWS"), which(mc == "CXRS"), which(mc == "CXSWS"))
indexEra <- c(which(era == "Pre-1920"), which(era == "1921-1940"), which(era == "1941-1960"), 
              which(era == "1961-1980"), which(era == "1981-2000"), which(era == "2001-2016"))
indexBP <- c(which(bp == "AAFC BEAVERLODGE RS"), which(bp == "AAFC CEREAL RC/WINNIPEG"), which(bp == "AAFC LACOMBE RC"), 
             which(bp == "AAFC LETHBRIDGE RC"), which(bp == "AGRIPRO/SYNGENTA"), which(bp == "CROP DEVELOPMENT CENTRE"),
             which(bp == "SASKATCHEWAN WHEAT POOL"), which(bp == "SPA RC/AAFC SWIFT CURRENT"), which(bp == "UNIVERSTIY OF ALBERTA"),
             which(bp == "UNIVERSITY OF  MANITOBA"))
indexDesig <- c(which(desig == "HRS"), which(desig == "HRW"), which(desig == "HWS"), which(desig == "HWW"), 
                which(desig == "SRS"), which(desig == "SRW"), which(desig == "SWS"), which(desig == "SWW"))
indexHard <- c(which(desig == "HRS"), which(desig == "HWS"), which(desig == "HRW"), which(desig == "HWW")) # HRS HRW
indexSoft <- c(which(desig == "SRS"), which(desig == "SWS"), which(desig == "SRW"), which(desig == "SWW")) # SWS
indexRed <- c(which(desig == "HRS"), which(desig == "SRS"), which(desig == "HRW"), which(desig == "SRW")) # HRS HRW
indexWhite <- c(which(desig == "HWS"), which(desig == "SWS"), which(desig == "HWW"), which(desig == "SWW")) # SWS
indexSpring <- c(which(desig == "HRS"), which(desig == "HWS"), which(desig == "SRS"), which(desig == "SWS")) # HRS SWS
indexWinter <- c(which(desig == "HRW"), which(desig == "HWW"), which(desig == "SRW"), which(desig == "SWW")) # HRW

table(desig)

# load("Data\\Intermediate\\CCP\\full_kmeans.RDATA")
# CCP6 <- c("Group 6", "Group 5", "Group 1", "Group 4", "Group 3", "Group 2")[bests]
# CCP2 <- c("Group 2", "Group 2", "Group 1", "Group 2", "Group 1", "Group 1")[bests]

load("Data\\Intermediate\\dbscan\\wheatHdbscan.RData")
# wheatHdbscan$cluster[which(wheatHdbscan$cluster == 0)] <- wheatHdbscan$cluster[which(wheatHdbscan$cluster == 0)] + 5
wheatHdbscan6 <- factor(wheatHdbscan$cluster)
wheatHdbscan2 <- as.integer(wheatHdbscan6)
wheatHdbscan2[which(wheatHdbscan2 != 6)] <- 1

era2 <- as.character(era)
era2[-indexEra] <- "NA"
bp2 <- as.character(bp)
bp2[-indexBP] <- "NA"
desig2 <- as.character(desig)
desig2[-indexDesig] <- "NA"

desigSW <- as.character(desig)
desigSW[indexWinter] <- "winter"
desigSW[indexSpring] <- "spring"
desigSW[-c(indexWinter, indexSpring)] <- "NA"
desigRW <- as.character(desig)
desigRW[indexRed] <- "red"
desigRW[indexWhite] <- "white"
desigRW[-c(indexRed, indexWhite)] <- "NA"
desigHS <- as.character(desig)
desigHS[indexHard] <- "hard"
desigHS[indexSoft] <- "soft"
desigHS[-c(indexHard, indexSoft)] <- "NA"
desigSWRW <- paste(desigSW, desigRW)
desigSWHS <- paste(desigSW, desigHS)
desigRWHS <- paste(desigSW, desigRW)
desigHRSOther <- as.character(desig)
desigHRSOther[which(desigHRSOther != "HRS")] <- "other"


strata <- data.frame("SW" = as.factor(desigSW), "RW" = as.factor(desigRW), "HS" = as.factor(desigHS),
                     "SWRW" = as.factor(desigSWRW), "SWHS" = as.factor(desigSWHS), "RWHS" = as.factor(desigRWHS),
                     "HRSOther" = as.factor(desigHRSOther), "desig" = as.factor(desig2), "era" = as.factor(era2), "bp"= as.factor(bp2),
                     "hdbscan2" = wheatHdbscan2, "hdbscan6" = wheatHdbscan6, "dend" = clusters)

genind <- df2genind(t(data.frame(genotypes)),
                    ind.names = as.character(sample.id),
                    loc.names = snp.id, NA.char = "N", ploidy = 1, type = "codom",
                    ncode = 1, strata = strata)

save(genind, file = paste0("Data\\Intermediate\\Adegenet\\genind_wheat.RData"))

# ## loading the gds of the data and pullling some attributes out
# gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds"
# source("src\\functions\\data_loading.R")
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
