library(SNPRelate)
library(ape)
library(circlize)
library(dendextend)
library(extrafont)
library(RColorBrewer)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the population data
wheat <- snpgdsOpen("Data\\Formatted\\HRS_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.chr <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
meta.data <- cbind(BP = as.character(read.gdsn(index.gdsn(wheat, "samp.annot/BP"))),
                   Year = as.character(read.gdsn(index.gdsn(wheat, "samp.annot/Year"))),
                   MC = as.character(read.gdsn(index.gdsn(wheat, "samp.annot/MC"))),
                   Desig = as.character(read.gdsn(index.gdsn(wheat, "samp.annot/designation"))))
rownames(meta.data) <- sample.id

# reformat some stuff
desig <- split(sample.id, meta.data[,3])

# neighbour joining
neighbour <- nj(1 - snpgdsIBS(wheat, autosome.only = F)$ibs)
labels(neighbour) <- sample.id
neighbour <- root(neighbour, which(sample.id == "Selkirk"), r = TRUE)
neighbour <- chronos(neighbour, lambda = 0, model = "correlated")
neighbour <- as.dendrogram(neighbour)

## printing out the neighbourrogram
png("Analysis\\Figures\\dend\\HRS_nj_dend.png", 
    family="Times New Roman", width = 11, height = 11, pointsize = 20, units = "in", res = 500)

circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, track.margin = c(0.005, 0.005))
circos.initialize("foo", xlim = c(0, 234), sector.width = 1)

##
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:234, rep(0, 234), labels(neighbour),
              facing = "clockwise", niceFacing = T, cex = 0.375, adj = c(0, -0.5), font = 2)
}, track.height = 0.15, bg.border = NA)

##
pop.code <- cut(as.numeric(as.character(read.gdsn(index.gdsn(wheat, "samp.annot/Year")))),
                breaks = c(1800, 1900, 1920, 1940, 1960, 1980, 2000, 2020), labels = F)
pop.code[is.na(pop.code)] <- 8
pop.code <- factor(pop.code)
colours <- colors()[c(554, 53, 652, 48, 109, 26, 547, 153)]
legend("topright", legend=c("Pre-1900", "1901-1920", "1921-1940", "1941-1960", "1961-1980", "1981-2000",
                            "2001-2016", "Unknown"), pch=19, col=colours, cex = 0.45, bg = colors()[468],
       title = "Period of Release")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:233.05, rep(0,234), 0.95:233.95, rep(1,234), border = colors()[468], lwd = 0.8,
              col = sapply(as.numeric(as.character(pop.code))[match(labels(neighbour), sample.id)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)

##
pop.code <- read.gdsn(index.gdsn(wheat, "samp.annot/BP"))
colours <- colors()[c(116, 554, 53, 26, 100, 490, 42, 115, 139, 153, 451, 652, 114, 195, 557, 79, 48)]
legend("bottomright", legend=levels(pop.code), pch=19, col=colours, cex = 0.35, bg = colors()[142],
       title = "Breeding Program")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:233.05, rep(0,234), 0.95:233.95, rep(1,234), border = colors()[142], lwd = 0.8,
              col = sapply(as.numeric(pop.code)[match(labels(neighbour), sample.id)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)

##
pop.code <- read.gdsn(index.gdsn(wheat, "samp.annot/MC"))
colours <- colors()[c(195, 26, 79, 652, 419, 554, 153)]
legend("bottomleft", legend=levels(pop.code), pch=19, col=colours, cex = 0.45, bg = colors()[105],
       title = "Market Class")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:233.05, rep(0,234), 0.95:233.95, rep(1,234), border = colors()[105], lwd = 0.8,
              col = sapply(as.numeric(pop.code)[match(labels(neighbour), sample.id)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)

##
load("Data\\Formatted\\hrs_kmeans.RDATA")
pop.code <- bests[match(labels(neighbour), sample.id)]
colours <-  colors()[c(554, 100, 652, 53)]
# pie(rep(1,4), col = colours)
legend("topleft", legend=c("Group 1", "Group 2", "Group 3", "Group 4"),
       pch=19, col=colors()[c(554, 100, 53, 652)], cex = 0.4, bg = colors()[109], 
       title = "K-means Groups")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:233.05, rep(0,234), 0.95:233.95, rep(1,234), border = colors()[109], lwd = 0.8,
              col = colours[pop.code])
}, track.height = 0.04, bg.border = NA)

##
max_height <- attr(neighbour, "height")
neighbour <- color_branches(neighbour, k = 4, col = c(brewer.pal(3, "Dark2"), "black"))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(neighbour, max_height = max_height)
}, track.height = 0.40, bg.border = NA)
legend("center", title = "NJ Groups", legend=c("Group A", "Group B", "Group C"),
       pch=19, col=brewer.pal(3, "Dark2"), cex = 0.4)
# title("NJ Dendrogram with K-means Groups of HRS Varieties 
                                                    # with metadata", out = T, cex = 1, line = -2)
circos.clear()
dev.off()
snpgdsClose(wheat)
