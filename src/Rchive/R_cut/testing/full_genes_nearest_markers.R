library(SNPRelate)
library(plyr)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA")

wheat <-
  snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
snp.chrom <-
  as.factor(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
desig <- read.gdsn(index.gdsn(wheat, "samp.annot/designation"))
snpgdsClose(wheat)

load("Data\\Intermediate\\Aligned_genes\\wheat_genes_contigs.RData")

## finds the markers that are before and after the aligned gene's location if they exist
genes_markers <-
  by(cbind(snp.pos, genes_contigs, genotypes), snp.chrom, function (chr) {
    size <- dim(chr)[2]
    bf_afs <- list()
    for (i in 1:size) {
      if (!is.na(chr[i, 4])) {
        index <- which.min(abs(chr[, 1] - chr[i, 3]))
        if (chr[index, 1] > chr[i, 3]) {
          bf_af <- rbind(chr[index - 1, 5:size], chr[index, 5:size])
          if (dim(bf_af)[1] < 2) {
            row.names(bf_af) <- paste(chr[i, 4], index, sep = "_")
          } else {
            row.names(bf_af) <-
              c(paste(chr[i, 4], index - 1, sep = "_"),
                paste(chr[i, 4], index, sep = "_"))
          }
          bf_afs[[i]] <- bf_af
        } else {
          bf_af <- rbind(chr[index, 5:size], chr[index + 1, 5:size])
          row.names(bf_af) <-
            c(paste(chr[i, 4], index, sep = "_"),
              paste(chr[i, 4], index + 1, sep = "_"))
          bf_afs[[i]] <- bf_af
        }
      }
    }
    return(bf_afs)
  })
lapply(unlist(genes_markers, recursive = F), function (gene_markers) {
  row.names(gene_markers)[1]
})
index.HRS <- which(desig == "HRS")
index.HRW <- which(desig == "HRW")
index.SWS <- which(desig == "SWS")


## for the VRN homeologs identifies whether each is the same or different in the individuals of index1 compared to index2
vrn_homeos <- list()
count <- 1
# for (indexes in list(list(index.HRS, index.HRS), list(index.HRW, index.HRW), list(index.HRS, index.HRW))) {
for (indexes in list(list(index.HRS, index.HRS),
                     list(index.SWS, index.SWS),
                     list(index.HRS, index.SWS))) {
  vrn_homeo <- list()
  # for (homeo in c ("VRN-A1", "VRN-B1", "VRN-D1")) {
  for (homeo in c ("Tamyb10-A1", "Tamyb10-B1", "Tamyb10-D1")) {
    lapply(unlist(genes_markers, recursive = F), function (gene_markers) {
      tests <- list()
      if (grepl(homeo, row.names(gene_markers)[1])) {
        for (i in 1:length(indexes[[1]])) {
          if (3 %in% gene_markers[, i]) {
            next
          }
          tests[[i]] <- vector()
          for (j in 1:length(indexes[[2]])) {
            if (3 %in% gene_markers[, j]) {
              next
            }
            if (identical(gene_markers[, i], gene_markers[, j])) {
              tests[[i]][j] <- TRUE
            } else {
              tests[[i]][j] <- FALSE
            }
          }
        }
      }
      if (length(tests) != 0) {
        vrn_homeo[[homeo]] <<- tests
      }
    })
  }
  vrn_homeos[[count]] <- vrn_homeo
  count <- count + 1
}
# str(vrn_homeos)

## HRS, HRW, HRSvsHRW
lapply(vrn_homeos, function (vrn_homeo) {
  # str(vrn_homeo)
  dif <- 0
  same <- 0
  for (i in 1:length(vrn_homeo$`VRN-A1`)) {
    A <- vrn_homeo$`VRN-A1`[[i]]
    B <- vrn_homeo$`VRN-B1`[[i]]
    D <- vrn_homeo$`VRN-D1`[[i]]
    if (is.null(A) | is.null(B) | is.null(D)) {
      next
    }
    for (j in 1:length(A)) {
      if (is.na(A[j]) | is.na(B[j]) | is.na(D[j])) {
        next
      } else if (!A[j] | !B[j] | !D[j]) {
        dif <- dif + 1
      } else {
        same <- same + 1
      }
    }
  }
  
  return(list(dif / (dif + same), same / (dif + same)))
})

## HRS, HRW, HRSvsHRW
lapply(vrn_homeos, function (vrn_homeo) {
  # str(vrn_homeo)
  dif <- 0
  same <- 0
  for (i in 1:length(vrn_homeo$`Tamyb10-A1`)) {
    A <- vrn_homeo$`Tamyb10-A1`[[i]]
    B <- vrn_homeo$`Tamyb10-B1`[[i]]
    D <- vrn_homeo$`Tamyb10-D1`[[i]]
    if (is.null(A) | is.null(B) | is.null(D)) {
      next
    }
    for (j in 1:length(A)) {
      if (is.na(A[j]) | is.na(B[j]) | is.na(D[j])) {
        next
      } else if (!A[j] | !B[j] | !D[j]) {
        dif <- dif + 1
      } else {
        same <- same + 1
      }
    }
  }
  return(list(dif / (dif + same), same / (dif + same)))
})
