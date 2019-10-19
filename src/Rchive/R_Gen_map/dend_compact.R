library(SNPRelate)
library(circlize)
library(dendextend)
library(extrafont)
library(scrime)


## loading the data
wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset_imputed.gds")
# wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset.gds")
sample.id <- read.gdsn(index.gdsn(wheat.subset, "sample.id"))
genotypes <- read.gdsn(index.gdsn(wheat.subset, "genotype"))

## producing the dendrogram
wheat.IBS <- snpgdsIBS(wheat.subset)
dend <- as.dendrogram(hclust(as.dist(1-wheat.IBS$ibs), method = "complete"))

## printing out the dendrogram
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\dend\\dend_compact_complete.png", 
    family="Times New Roman", width = 11, height = 11, pointsize = 20, units = "in", res = 300)
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, track.margin = c(0.005, 0.005))
circos.initialize("foo", xlim = c(0, 359), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:359, rep(0, 359), as.character(sample.id)[labels(dend)],
              facing = "clockwise", niceFacing = T, cex = 0.25, adj = c(0, -0.1), font = 2)
}, track.height = 0.15, bg.border = NA)
##
pop.code <- cut(as.numeric(as.character(read.gdsn(index.gdsn(wheat.subset, "samp.annot/Year")))),
                          breaks = c(1800, 1900, 1920, 1940, 1960, 1980, 2000, 2020), labels = F)
pop.code[is.na(pop.code)] <- 8
pop.code <- factor(pop.code)
colours <- colors()[c(554, 53, 652, 48, 109, 26, 547, 153)]
legend(0.88, 0.8, legend=c("Pre-1900", "1901-1920", "1921-1940", "1941-1960", "1961-1980", "1981-2000", 
                            "2001-2016", "Unknown"), pch=19, col=colours, cex = 0.45, bg = colors()[468],
       title = "Period of Release")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:358.05, rep(0,359), 0.95:358.95, rep(1,359), border = colors()[468], lwd = 0.6,
              col = sapply(as.numeric(as.character(pop.code))[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)
##
pop.code <- read.gdsn(index.gdsn(wheat.subset, "samp.annot/BP"))
colours <- colors()[c(116, 554, 53, 26, 100, 12, 490, 42, 115, 139, 153, 451, 652, 114, 195, 557, 79, 48)]
legend("bottomright", legend=levels(pop.code), pch=19, col=colours, cex = 0.35, bg = colors()[142],
       title = "Breeding Program")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:358.05, rep(0,359), 0.95:358.95, rep(1,359), border = colors()[142], lwd = 0.6,
              col = sapply(as.numeric(pop.code)[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)
##
pop.code <- read.gdsn(index.gdsn(wheat.subset, "samp.annot/MC"))
colours <- colors()[c(195, 26, 49, 79, 652, 100, 419, 115, 114, 53, 554, 450, 153)]
legend("bottomleft", legend=levels(pop.code), pch=19, col=colours, cex = 0.45, bg = colors()[105],
       title = "Market Class")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:358.05, rep(0,359), 0.95:358.95, rep(1,359), border = colors()[105], lwd = 0.6,
              col = sapply(as.numeric(pop.code)[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)
##
pop.code <- read.gdsn(index.gdsn(wheat.subset, "samp.annot/designation"))
colours <- colors()[c(554, 115, 1, 26, 153, 450, 114, 652, 153)]
legend(-1.08, 0.8, legend=levels(pop.code), pch=19, col=colours, cex = 0.45, bg = colors()[525],
       title = "Designation")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:358.05, rep(0,359), 0.95:358.95, rep(1,359), border = colors()[525], lwd = 0.6,
              col = sapply(as.numeric(pop.code)[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)
##
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
palette(colours)
dend_colour <- color_branches(dend, k = 7, col = colours)
max_height <- attr(dend_colour, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(dend_colour, max_height = max_height)
}, track.height = 0.40, bg.border = NA)
legend("center", legend=c("Hard Red Spring", "Chinese", "Australian", 
                          "Prairie Spring Red", "White Spring", 
                          "Winter - Out", "Winter"), 
       pch=19, col=1:7, cex = 0.4)
title("Fig 2: Dendrogram of Wheat Varieties with Group Data", out = T, cex = 1.5, line = -1)
circos.clear()
dev.off()
snpgdsClose(wheat.subset)
