library(SNPRelate)
library(extrafont)
library(plyr) #using revalue function
library(data.table)
library(devEMF)

## Setting up
# Loading the data
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data_2.gds")
snp.id <- read.gdsn(index.gdsn(wheat.data, "snp.id"))
snp.pos <- read.gdsn(index.gdsn(wheat.data, "snp.position"))/1000
snp.chr <- factor(read.gdsn(index.gdsn(wheat.data, "snp.chromosome")), levels = 0:22)

# Removing SNPs with high MR and low MAF
snpset <- snpgdsSelectSNP(wheat.data, autosome.only = F, maf = 0.05, missing.rate = 0.1)
snpset.id <- unlist(snpset)
informative <- match(snpset.id, snp.id)
snp.pos.inf <- snp.pos[informative]
snp.chr.inf <- snp.chr[informative]
snp.id.inf <- snp.id[informative]
snpgdsClose(wheat.data)

# creating vectors without unplaced SNPs
snp.pos.inf.non.na <- snp.pos.inf[-which(snp.chr.inf == 0)]
snp.chr.inf.non.na <- droplevels(snp.chr.inf[-which(snp.chr.inf == 0)])
snp.id.inf.non.na <- snp.id.inf[-which(snp.chr.inf == 0)]
# re-order snp.chr factor levels
snp.chr.inf.non.na.reordered <- factor(snp.chr.inf.non.na, levels = c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3)))

# Chromosome labels
labels <- c("NA", "1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", 
            "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D")
# genome chromosome label order
labels.reordered <- labels[c(1, seq(2, 20, 3), seq(3, 21, 3), seq(4, 22, 3))]


#######################################################################################################
#######################################################################################################
## Plotting positon by index number
# some values that are reused on the graphing function
num.snps <- table(snp.chr.inf.non.na.reordered)
sum(num.snps[8:14])/sum(num.snps)
cex <- 0.8
xlim <- c(0, round_any(max(num.snps), 10))
ylim <- c(0, round_any(max(snp.pos.inf.non.na), 10, f = ceiling))
xlabs <- c(round_any(seq(0, max(num.snps), max(num.snps)/4), 10))
ylabs <- c(round_any(seq(0, max(snp.pos.inf.non.na), max(snp.pos.inf.non.na)/4), 10, f = ceiling))
# making the graph itself
pdf("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\snp_pos_x_index.pdf", family="Impact")
# making 21 graphs
par(mfrow = c(3,7), oma = c(5,5,2,1))
# counts how many times tapply is called
count = 1
# allows us to operate on markers from each chromosome separately
foo <- tapply(snp.pos.inf.non.na, snp.chr.inf.non.na.reordered, function (x) {
  par(mar = c(2, 0, 0, 0))
  # the first graph in each row
  if (count %% 7 == 1) {
    plot(x, pch = "-", ylim = ylim, xlim = xlim,
         las = 3, yaxt = "n", xaxt = "n", cex.axis = cex)
    axis(2, at = ylabs, labels = ylabs)
    # the first two rows
    if (count < 15) mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
                          line = 0.5, cex = cex)
    # the last row
    if (count == 15 ) {
      par(mar = c(2, 0, 0, 0))
      axis(1, at = xlabs, labels = xlabs, las = 2)
      mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
            line = 3.1, cex = cex)
    }
  }
  # the last 6 graphs in each row
  if (count %% 7 != 1) {
    plot(x, pch = "-", xlab = labels.reordered[count + 1], ylim = ylim, xlim = xlim,
         las = 3, yaxt = "n", xaxt = "n")
    # the first two rows
    if (count < 16) mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
                          line = 0.5, cex = cex)
    # the last row
    if (count > 15) {
      par(mar = c(2, 0, 0, 0))
      axis(1, at = xlabs, labels = rep("" ,5))
      mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
            line = 1.2, cex = cex)
    }
  }
  count <<- count + 1
})
# titles and axes labels
title("21,603 SNPs Ordered and Plotted by Postion on each Chromosome", out = T)
title(xlab = "Ordered Index Number of Markers", ylab = "Distance from Start of Chrosmosome in cM", 
      outer = T, cex.lab = 1.5)
dev.off()

## distance between adjacent positions
snp.diffs <- data.frame()
count <- 0
blah <- tapply(snp.pos.inf.non.na, snp.chr.inf.non.na.reordered, function (x) {
  count <<- count + 1
  snp.diffs <<- rbind(snp.diffs, cbind(count, diff(x)))
})
snp.diffs[,1] <- factor(snp.diffs[,1])

# counts the number of times each distance postion shows up, lets us calcluate the percent that are 0
table <- table(snp.diffs[,2])
sum(table[2:length(table)])/table[1]

#plots those distances
pdf("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\snp_distances_between_adjacent_hist.pdf",
    family="Impact")
plot(table, main = "Histogram of Distances between Adjacent Positions", 
     xlab = "Distance Between Unique Markers in cM", ylab = "Count")
dev.off()

# creates a graph of boxplots for the distance distribution of each chromosome
pdf("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\snp_distances_between_adjacent.pdf", family="Impact")
boxplot(snp.diffs[,2] ~ snp.diffs[,1], outline = F, 
        main = "Boxplots of Distance Distributions Between Postions on Chromsomes", 
        xlab = "Chromosome", ylab = "Distance Between Adjacent Markers in cM", xaxt = "n")
axis(1, at = 1:21, labels = labels.reordered[2:22], las = 2)
dev.off()

# finds the summary of distance distributions for each chromosome and turns them into a nice table
summary.by.chr <- data.frame()
blah <- tapply(snp.diffs[,2], snp.diffs[,1], function (x) {
  summary.by.chr <<- rbind(summary.by.chr, as.vector(summary(x)))
})
summary.by.chr <- rbind(summary.by.chr, as.vector(summary(snp.diffs[,2])))
row.names(summary.by.chr) <- c(labels.reordered[2:22], "All")
colnames(summary.by.chr) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
write.csv(summary.by.chr,
          "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\summary_snp_distances_between_adjacent.csv")


## Percent unique locations on chromosomes
length(unlist(tapply(snp.pos.inf.non.na, snp.chr.inf.non.na, unique)))/length(snp.pos.inf.non.na)

#######################################################################################################
#######################################################################################################
## graphing unique positions on chromosomes
# finds indices of unique positions on each chromosome
unique.indices <- data.frame()
length <- 0
blah <- tapply(snp.pos.inf.non.na, snp.chr.inf.non.na, function (x) {
  unique.indices <<- c(unique.indices, match(unique(x), x) + length)
  length <<- length + length(x)
})
unique.indices <- as.vector(unlist(unique.indices))
num.snps <- table(snp.chr.inf.non.na.reordered[unique.indices])
# some values that are reused on the grpahing function
cex <- 0.8
xlim <- c(0, round_any(max(num.snps), 10))
ylim <- c(0, round_any(max(snp.pos.inf.non.na), 5, f = ceiling))
xlabs <- c(round_any(seq(0, max(num.snps), max(num.snps)/4), 10))
ylabs <- c(round_any(seq(0, max(snp.pos.inf.non.na), max(snp.pos.inf.non.na)/4), 10, f = ceiling))
# making the graph itself
pdf("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\snp_pos_unique_x_index.pdf", 
    family="Impact")
# making 21 graphs
par(mfrow = c(3,7), oma = c(5,5,2,1))
# counts how many times tapply is called
count = 1
# allows us to operate on markers form each chromosome separately
blah <- tapply(snp.pos.inf.non.na[unique.indices], snp.chr.inf.non.na.reordered[unique.indices], 
               function (x) {
  par(mar = c(2, 0, 0, 0))
  # the first graph in each row
  if (count %% 7 == 1) {
    plot(x, pch = "--", ylim = ylim, xlim = xlim,
         las = 3, yaxt = "n", xaxt = "n", cex.axis = cex)
    axis(2, at = ylabs, labels = ylabs)
    # the first two rows
    if (count < 15) mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
                          line = 0.5, cex = cex)
    # the last row
    if (count == 15 ) {
      par(mar = c(2, 0, 0, 0))
      axis(1, at = xlabs, labels = xlabs, las = 2)
      mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
            line = 3.1, cex = cex)
    }
  }
  # the last 6 graphs in each row
  if (count %% 7 != 1) {
    plot(x, pch = "-", xlab = labels.reordered[count + 1], ylim = ylim, xlim = xlim,
         las = 3, yaxt = "n", xaxt = "n")
    # the first two rows
    if (count < 16) mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
                          line = 0.5, cex = cex)
    # the last row
    if (count > 15) {
      par(mar = c(2, 0, 0, 0))
      axis(1, at = xlabs, labels = rep("" ,5))
      mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
            line = 1.2, cex = cex)
    }
  }
  count <<- count + 1
})
# titles and axes labels
# title(paste(sum(num.snps), "Unique Loci Ordered and Plotted by Postion on Each Chromosome"), out = T)
title(xlab = "Ordered Index Number of Locus", ylab = "Distance from Start of Chrosmosome in cM",  
      outer = T, cex.lab = 1.5)
dev.off()


## distance between adjacent unique positions
snp.diffs <- data.frame()
count <- 0
blah <- tapply(snp.pos.inf.non.na[unique.indices], snp.chr.inf.non.na.reordered[unique.indices], 
               function (x) {
  count <<- count + 1
  snp.diffs <<- rbind(snp.diffs, cbind(count, diff(x)))
})
snp.diffs[,1] <- factor(snp.diffs[,1])

# counts the number of times each distance postion shows up
table <- table(round(snp.diffs[,2], 2))

#plots those distances
pdf("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\snp_distances_between_unique_hist.pdf",
    family="Impact")
par(mar = c(4,4,0,2), oma = c(0,0,2,0))
plot(table, xlab = "Distance Between Adjacent Unique Loci in cM", ylab = "Count")
dev.off()

# creates a graph of boxplots for the distance distribution of each chromosome
pdf("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\snp_distances_between_unique.pdf", 
    family="Impact")
par(mar = c(4,4,0,2), oma = c(0,0,2,0))
boxplot(snp.diffs[,2] ~ snp.diffs[,1], outline = F, 
        xlab = "Chromosome", ylab = "Distance Between Unique Adjacent Loci in cM", xaxt = "n")
axis(1, at = 1:21, labels = labels.reordered[2:22], las = 2)
dev.off()

# finds the summary of distance distributions for each chromosome and turns them into a nice table
summary.by.chr <- data.frame()
blah <- tapply(snp.diffs[,2], snp.diffs[,1], function (x) {
  summary.by.chr <<- rbind(summary.by.chr, as.vector(summary(x)))
})
summary.by.chr <- rbind(summary.by.chr, as.vector(summary(snp.diffs[,2])))
row.names(summary.by.chr) <- c(labels.reordered[2:22], "All")
colnames(summary.by.chr) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
write.csv(summary.by.chr, 
"C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\summary_snp_distances_between_unique.csv")

# summary of median values
summary(summary.by.chr[1:21,3])


#######################################################################################################
#######################################################################################################
## Looking at pruned SNP distances and distributions
# finds indices of unique positions on each chromosome
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data_2.gds")
set.seed(1000)
snp.set <- snpgdsLDpruning(wheat.data, ld.threshold=0.3, autosome.only = F,
                           maf = 0.05, missing.rate = 0.1, slide.max.bp = 0.345)
snpgdsClose(wheat.data)
snp.set.id <- unlist(snp.set)
kept.indices.non.na <- match(snp.set.id, snp.id.inf.non.na)
kept.indices.non.na <- kept.indices.non.na[1:4075]
num.snps <- table(snp.chr.inf.non.na.reordered[kept.indices.non.na])
# some values that are reused on the grpahing function
cex <- 0.8
xlim <- c(0, round_any(max(num.snps), 10))
ylim <- c(0, round_any(max(snp.pos.inf.non.na), 5, f = ceiling))
xlabs <- c(round_any(seq(0, max(num.snps), max(num.snps)/4), 10))
ylabs <- c(round_any(seq(0, max(snp.pos.inf.non.na), max(snp.pos.inf.non.na)/4), 10, f = ceiling))
# making the graph itself
pdf("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\snp_pos_kept_x_index.pdf")
# making 21 graphs
par(mfrow = c(3,7), oma = c(5,5,2,1))
# counts how many times tapply is called
count = 1
# allows us to operate on markers form each chromosome separately
blah <- tapply(snp.pos.inf.non.na[kept.indices.non.na], snp.chr.inf.non.na.reordered[kept.indices.non.na], function (x) {
  par(mar = c(2, 0, 0, 0))
  # the first graph in each row
  if (count %% 7 == 1) {
    plot(unique(x), pch = "-", ylim = ylim, xlim = xlim,
         las = 3, yaxt = "n", xaxt = "n", cex.axis = cex)
    axis(2, at = ylabs, labels = ylabs)
    # the first two rows
    if (count < 15) mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
                          line = 0.5, cex = cex)
    # the last row
    if (count == 15 ) {
      par(mar = c(2, 0, 0, 0))
      axis(1, at = xlabs, labels = xlabs, las = 2)
      mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
            line = 3.1, cex = cex)
    }
  }
  # the last 6 graphs in each row
  if (count %% 7 != 1) {
    plot(unique(x), pch = "-", xlab = labels.reordered[count + 1], ylim = ylim, xlim = xlim,
         las = 3, yaxt = "n", xaxt = "n")
    # the first two rows
    if (count < 16) mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
                          line = 0.5, cex = cex)
    # the last row
    if (count > 15) {
      par(mar = c(2, 0, 0, 0))
      axis(1, at = xlabs, labels = rep("" ,5))
      mtext(1, text = paste(labels.reordered[count + 1], ": ", num.snps[count]),
            line = 1.2, cex = cex)
    }
  }
  count <<- count + 1
})
# titles and axes labels
# title(paste(sum(num.snps),"Kept SNPs Ordered and Plotted by Postion on each Chromosome"), out = T)
title(xlab = "Ordered Index Number of Markers", ylab = "Distance from Start of Chrosmosome in cM",  
      outer = T, cex.lab = 1.5)
dev.off()

## distance between adjacent unique positions
snp.diffs <- data.frame()
count <- 0
blah <- tapply(snp.pos.inf.non.na[kept.indices.non.na], snp.chr.inf.non.na.reordered[kept.indices.non.na],
               function (x) {
  count <<- count + 1
  snp.diffs <<- rbind(snp.diffs, cbind(count, diff(x)))
})
snp.diffs[,1] <- factor(snp.diffs[,1])

# counts the number of times each distance postion shows up
table <- table(round(snp.diffs[,2], 2))

#plots those distances
pdf("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\snp_distances_between_kept_hist.pdf")
# zeroing 0 count
table[1] <- 0
par(mar=c(5, 4, 1, 2))
plot(table, xlab = "Distance Between Adjacent Kept Loci in cM", ylab = "Count")
dev.off()

# creates a graph of boxplots for the distance distribution of each chromosome
pdf("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\snp_distances_between_kept.pdf")
par(mar=c(5, 4, 1, 2))
boxplot(snp.diffs[,2] ~ snp.diffs[,1], outline = F, 
        # main = "Boxplots of Distance Distributions Between Kept Postions on Chromsomes", 
        xlab = "Chromosome", ylab = "Distance Between Adjacent Kept SNPs in cM", xaxt = "n")
axis(1, at = 1:21, labels = labels.reordered[2:22], las = 2)
dev.off()

# finds the summary of distance distributions for each chromosome and turns them into a nice table
summary.by.chr <- data.frame()
blah <- tapply(snp.diffs[,2], snp.diffs[,1], function (x) {
  summary.by.chr <<- rbind(summary.by.chr, as.vector(summary(x)))
})
summary.by.chr <- rbind(summary.by.chr, as.vector(summary(snp.diffs[,2])))
row.names(summary.by.chr) <- c(labels.reordered[2:22], "All")
colnames(summary.by.chr) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
write.csv(summary.by.chr, "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\distances\\summary_by_chr_kept.csv")

# summary of median values
summary(summary.by.chr[1:21,3])