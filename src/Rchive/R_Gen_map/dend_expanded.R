library(SNPRelate)
library(circlize)
library(dendextend)
library(extrafont)
library(RColorBrewer)
library(Hmisc)
library(pvclust)
library(devEMF)

## loading the pvclust cluster object
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes_pruned_pvclust.RDATA")
dend <- genotypes.pvclust %>% as.dendrogram

## loading the population data
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data_2.gds")
sample.id <- read.gdsn(index.gdsn(wheat.data, "sample.id"))

## printing out the dendrogram
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\dend\\dend_expanded.png", 
    family="Times New Roman", width = 12, height = 12, pointsize = 20, units = "in", res = 300)
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, track.margin = c(0.005, 0.005))
circos.initialize("foo", xlim = c(0, 388), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:388, rep(0, 388), as.character(sample.id)[labels(dend)],
              facing = "clockwise", niceFacing = T, cex = 0.25, adj = c(0, -0.1), font = 2)
}, track.height = 0.15, bg.border = NA)
pop.code <- read.gdsn(index.gdsn(wheat.data, "origin"))
colours <- colors()[c(450, 139, 419, 554, 144, 53, 115, 153, 26)]
legend("topright", legend=levels(pop.code), pch=19, col=colours, cex = 0.45, bg = colors()[468], 
       title = "Country of Breeding")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:387.05, rep(0,388), 0.95:387.95, rep(1,388), border = colors()[468], lwd = 0.6,
              col = sapply(as.numeric(pop.code)[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)
pop.code <- read.gdsn(index.gdsn(wheat.data, "strength"))
colours <- colors()[c(376, 153, 424, 153)]
legend("bottomright", legend=levels(pop.code), pch=19, col=colours, cex = 0.40, bg = colors()[142],
       title = "Strength")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:387.05, rep(0,388), 0.95:387.95, rep(1,388), border = colors()[142], lwd = 0.6,
              col = sapply(as.numeric(pop.code)[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)
pop.code <- read.gdsn(index.gdsn(wheat.data, "colour"))
colours <- colors()[c(153, 554, 153, 1)]
legend("bottomleft", legend=levels(pop.code), pch=19, col=colours, cex = 0.45, bg = colors()[105],
       title = "Kernel Colour")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:387.05, rep(0,388), 0.95:387.95, rep(1,388), border = colors()[105], lwd = 0.6,
              col = sapply(as.numeric(pop.code)[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)
pop.code <- read.gdsn(index.gdsn(wheat.data, "season"))
colours <- colors()[c(153, 48, 153, 26)]
legend("topleft", legend=levels(pop.code), pch=19, col=colours, cex = 0.45, bg = colors()[525],
       title = "Planting Season")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:387.05, rep(0,388), 0.95:387.95, rep(1,388), border = colors()[525], lwd = 0.6,
              col = sapply(as.numeric(pop.code)[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
palette(colours)
dend_colour <- color_branches(dend, k = 7, col = colours)
max_height <- attr(dend_colour, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(dend_colour, max_height = max_height)
}, track.height = 0.40, bg.border = NA)
legend("center", legend=c("Red Spring", "Chinese", "Australian", 
                          "Extra Strong/\nPrairie", "White Spring", 
                          "Red Winter -\nOut", "Red Winter"), 
       pch=19, col=1:7, cex = 0.4)
text(0, 1, bquote(underline("Major Phenotype data of Wheats and Dendrogram Made by UPGMA")))
circos.clear()
dev.off()
snpgdsClose(wheat.data)