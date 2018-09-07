library(tidyverse)
library(SNPRelate)
library(adegenet)
library(dbscan)

## loading the gds of the data and pullling some attributes out
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\functions\\data_loading.R")

cluster <- read_rds("Data\\Intermediate\\dbscan\\full_hdbscan.rds")$cluster

# making indexes in order to create subsets containg just the varieties of the
# phenotype/cluster groups we are interested in
hrs <- which(desig == "HRS")
index_chrs <- hrs[which(hrs %in% which(cluster == 5))]
sws <- which(desig == "SWS")
index_csws <- sws[which(sws %in% which(cluster == 2))]
index_hrw <- which(desig == "HRW")

index_chrs_csws <- c(index_chrs, index_csws)
index_chrs_hrw <- c(index_chrs, index_hrw)
index_csws_hrw <- c(index_csws, index_hrw)

grouping <- list(list(index_chrs_csws, index_chrs_hrw, index_csws_hrw),
                 list("CHRS_CSWS", "CHRS_HRW", "HRW_CSWS"))

desig <- as.character(desig)
for (i in 1:length(grouping)) {
  index <- grouping[[1]][[i]]
  name <- grouping[[2]][[i]]
  
  strata <- data.frame(desig[index])
  colnames(strata) <- name
  
  genind <- df2genind(t(data.frame(genotypes[,index])),
                      ind.names = as.character(sample_id)[index],
                      loc.names = snp_id, NA.char = "N", ploidy = 1, 
                      type = "codom", ncode = 1, strata = strata)
  
  write_rds(genind, 
       path = paste0("Data\\Intermediate\\Adegenet\\", name, "_genind.rds"))
}

# making subsetted geninds of groups containing and not containing alleles of
# certain lR genes
gene_pres <- read_csv(
    "Data\\Intermediate\\Aligned_genes\\gene_presence_randhawa.csv", 
    col_names = c("sample", "gene_A", "gene_B", "gene_C"))
for (gene in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
  indivs <- c(which(gene_pres$gene_A == gene), which(gene_pres$gene_B == gene),
              which(gene_pres$gene_C == gene))
  index_indivs <- match(gene_pres$sample[indivs], sample_id)
  index_indivs <- index_indivs[!is.na(index_indivs)]
  index_not_indivs <- match(gene_pres$sample[-indivs], sample_id)
  index_not_indivs <- index_not_indivs[!is.na(index_not_indivs)]

  strata <- data.frame(c(rep(gene, length(index_indivs)),
                       rep(paste("not", gene), length(index_not_indivs))))
  colnames(strata) <- gene

  genind <- df2genind(t(data.frame(
                        genotypes[,c(index_indivs, index_not_indivs)])),
                      ind.names = as.character(
                                  sample_id)[c(index_indivs,
                                               index_not_indivs)],
                      loc.names = snp_id, NA.char = "N", ploidy = 1,
                      type = "codom", ncode = 1, strata = strata)

  write_rds(genind, 
       path = paste0("Data\\Intermediate\\Adegenet\\", gene, "_genind.rds"))
}