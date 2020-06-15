library(SNPRelate)
library(randomForest)
library(circlize)
library(extrafont)
library(gtools)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\R\\functions\\colour_sets.R")
source("Analysis\\R\\functions\\funcs_calc_IBD.R")
source("Analysis\\R\\functions\\funcs_draw_links.R")

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds"
source("Analysis\\R\\functions\\data_loading.R")

load("Data\\Intermediate\\CCP\\CCP.RDATA")
labelOrder <- results[[6]]$consensusTree$order

load("Data\\Intermediate\\CCP\\full_kmeans.RDATA")
bests <- as.factor(bests)
levels(bests) <- c("Group 6", "Group 5", "Group 1", "Group 4", "Group 2","Group 3")

# genotypes <- replace(genotypes, genotypes == 3, NA)
# genotypes <- replace(genotypes, genotypes == 0, 1)

# genotypes.imputed <- rfImpute(bests ~ ., genotypes)
# save(genotypes.imputed, file = "Data\\Intermediate\\genotypes_imputed.RData")
load("Data\\Intermediate\\genotypes_imputed.RData")

# totalDistances <- hapLengthsTotalFinder(windowMaker(snp.pos, genotypes.imputed, snp.chrom))
# save(totalDistances, file = "Data\\Intermediate\\haplotype_sharing\\total_distances.RData")
load("Data\\Intermediate\\haplotype_sharing\\total_distances.RData")

totalDist <- as.dist(t(totalDistances))

chromLengths <- by(snp.pos, snp.chrom, function (x) {
  x[length(x)] - x[1]
})
sumChromLengths <- sum(as.numeric(chromLengths))
percentIBD <- totalDistances/sumChromLengths*100

sq <- seq(0.9, 0.99, 0.01)

## links between HRS groups
apply(permutations(3, 2, c(1, 2, 4)), 1, function (i) {
  if (i[1] < i[2]) {
    createLinkPlot("HRS", sq, 0.01, i[1], i[2])
  }
})

## links within groups
apply(cbind(1:6, 1:6), 1, function (i) {
    createLinkPlot("self", sq, 0.01, i[1], i[2])
})

sq <- c(0.9, 0.95)

## links between non-hrs groups
apply(permutations(3, 2, c(3, 5, 6)), 1, function (i) {
  if (i[1] < i[2]) {
    createLinkPlot("non-HRS", sq, 0.05, i[1], i[2])
  }
})

## links between HRS and non-HRS groups
apply(expand.grid(c(1, 2, 4), c(3, 5, 6)), 1, function (i) {
    createLinkPlot("between", sq, 0.05, i[1], i[2])
})

## links between indivs
for (i in 1:358) {
  createLinkPlot("indiv", 0.90, 0.10, i, "all", plotDesig = TRUE)
}

## find linked indivs to each indiv
indivsIndivs <- list()
for (i in 1:358) {
  indivsIndivs[[i]] <- specificLinks(0.90, 0.10, i)
}
indivsIndivs[[331]]
byLinks <- sample.id[labelOrder][Reduce(intersect, indivsIndivs[which(counts >= quantile(counts, 0.97))])]

source("Analysis\\R\\functions\\funcs_draw_links.R")
## find average link freq of inidvs link in top x percent
indivsAves <- vector(length = 358)
for (i in 1:358) {
  indivsAves[[i]] <- averageLinks(0.90, 0.10, i)
}
min(indivsAves, na.rm = T)
which(indivsAves == max(indivsAves, na.rm = T))
byAve <- sample.id[labelOrder][Reduce(intersect, indivsIndivs[which(indivsAves >= quantile(indivsAves, 0.99, na.rm = T))])]

intersect(byLinks, byAve)

## linkage freq between groups in top 10%
apply(expand.grid(1:6, 1:6), 1, function (i) {
  if (i[1] >= i[2]) {
    linkage(0.90, 1, i[1], i[2])
  }
})


## finding average linakge stat
apply(expand.grid(1:6, 1:6), 1, function (group) {
  lengths <- 0
  count <- 0
  if (group[1] >= group[2]) {
    for (i in 1:ncol(totalDistances)) {
      for (j in i+1:ncol(totalDistances)) if (j < ncol(totalDistances)) {
        first <- as.numeric(bests[labelOrder])[i]
        second <- as.numeric(bests[labelOrder])[j]
        if (first == group[1] & second == group[2] |
            first == group[2] & second == group[1]) {
          lengths <- c(lengths, percentIBD[i,j])
          count <- count + 1
        }
      }
    }
  }
  paste("The average length of shared haplotypes between groups",
        levels(bests)[group[1]], "and", levels(bests)[group[2]],
        "is", sum(lengths)/count, 
        "its standard deviation is", sd(lengths))
})


## average IBD between indivs of groups 1, 2, & 3
lengths2 <- 0
count <- 0
for (i in 1:ncol(totalDistances)) {
  for (j in i+1:ncol(totalDistances)) if (j < ncol(totalDistances)) {
    first <- bests[labelOrder][i]
    second <- bests[labelOrder][j]
    if ((first == "Group 1" & second == "Group 2" | first == "Group 2" & second == "Group 1") |
        (first == "Group 1" & second == "Group 3" | first == "Group 3" & second == "Group 1") |
        (first == "Group 2" & second == "Group 3" | first == "Group 3" & second == "Group 2")) {
      lengths2 <- c(lengths2, percentIBD[i,j])
      count <- count + 1
    }
  }
}
print(paste("The average length of shared haplotypes between groups 1, 2, and 3 is",
            sum(lengths2)/count, "its standard deviation is", sd(lengths2)))

## average IBD of individuals within groups 1, 2, and 3
lengths3 <- 0
count <- 0
for (i in 1:ncol(totalDistances)) {
  for (j in i+1:ncol(totalDistances)) if (j < ncol(totalDistances)) {
    first <- bests[labelOrder][i]
    second <- bests[labelOrder][j]
    if (first == "Group 1" & second == "Group 1" |
        first == "Group 2" & second == "Group 2" |
        first == "Group 3" & second == "Group 3" ) {
      lengths3 <- c(lengths3, percentIBD[i,j])
      count <- count + 1
    }
  }
}
print(paste("The average length of shared haplotypes within groups 1, 2, and 3 is",
            sum(lengths3)/count, "its standard deviation is", sd(lengths3)))

## avaerge IBD between individuals of groups 4, 5, and 6
lengths4 <- 0
count <- 0
for (i in 1:ncol(totalDistances)) {
  for (j in i+1:ncol(totalDistances)) if (j < ncol(totalDistances)) {
    first <- bests[labelOrder][i]
    second <- bests[labelOrder][j]
    if ((first == "Group 4" & second == "Group 5" | first == "Group 5" & second == "Group 4") |
        (first == "Group 4" & second == "Group 6" | first == "Group 6" & second == "Group 4") |
        (first == "Group 5" & second == "Group 6" | first == "Group 6" & second == "Group 5")) {
      lengths4 <- c(lengths4, percentIBD[i,j])
      count <- count + 1
    }
  }
}
print(paste("The average length of shared haplotypes between groups 4, 5, and 6 is",
            sum(lengths4)/count, "its standard deviation is", sd(lengths4)))

t.test(lengths, lengths4)

# average IBD across groups 4, 5, and 6 of individuals within groups 4, 5, and 6
lengths5 <- 0
count <- 0
for (i in 1:ncol(totalDistances)) {
  for (j in i+1:ncol(totalDistances)) if (j < ncol(totalDistances)) {
    first <- bests[labelOrder][i]
    second <- bests[labelOrder][j]
    if (first == "Group 4" & second == "Group 4" |
        first == "Group 5" & second == "Group 5" |
        first == "Group 6" & second == "Group 6" ) {
      lengths5 <- c(lengths5, percentIBD[i,j])
      count <- count + 1
    }
  }
}
print(paste("The average length of shared haplotypes within groups 4, 5, and 6 is",
            sum(lengths5)/count, "its standard deviation is", sd(lengths5)))

t.test(lengths3, lengths5)

lengths6 <- 0
count <- 0
for (i in 1:ncol(totalDistances)) {
  for (j in i+1:ncol(totalDistances)) if (j < ncol(totalDistances)) {
    first <- bests[labelOrder][i]
    second <- bests[labelOrder][j]
    if ((first == "Group 1" & second == "Group 4" | first == "Group 4" & second == "Group 1") |
        (first == "Group 1" & second == "Group 5" | first == "Group 5" & second == "Group 1") |
        (first == "Group 1" & second == "Group 6" | first == "Group 6" & second == "Group 1") |
        (first == "Group 2" & second == "Group 4" | first == "Group 4" & second == "Group 2") |
        (first == "Group 2" & second == "Group 5" | first == "Group 5" & second == "Group 2") |
        (first == "Group 2" & second == "Group 6" | first == "Group 6" & second == "Group 2") |
        (first == "Group 3" & second == "Group 4" | first == "Group 4" & second == "Group 3") |
        (first == "Group 3" & second == "Group 5" | first == "Group 5" & second == "Group 3") |
        (first == "Group 3" & second == "Group 6" | first == "Group 6" & second == "Group 3")) {
      lengths6 <- c(lengths7, percentIBD[i,j])
      count <- count + 1
    }
  }
}
print(paste("The average length of shared haplotypes between groups 1, 2, & 3 and 4, 5, & 6 is",
            sum(lengths6)/count, "its standard deviation is", sd(lengths6)))




# red 1
# dark red 2
# orange 4
# yellow 5
# teal 3
# purple 6