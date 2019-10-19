library(SNPRelate)
library(circlize)
library(dendextend)
library(extrafont)
library(ape)
library(RColorBrewer)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the population data
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))

# neighbour joining
neighbour <- nj(1 - snpgdsIBS(wheat, autosome.only = F)$ibs)
labels(neighbour) <- sample.id
neighbour <- root(neighbour, which(sample.id == "Selkirk"), r = TRUE)
neighbour <- chronos(neighbour, lambda = 0, model = "correlated")
neighbour <- as.dendrogram(neighbour)

## printing out the neighbourrogram
png("Results\\dend\\full_nj_dend.png", 
    family="Times New Roman", width = 11, height = 11, pointsize = 20, units = "in", res = 500)

circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, track.margin = c(0.005, 0.005))
circos.initialize("foo", xlim = c(0, 358), sector.width = 1)

##
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:358, rep(0, 358), labels(neighbour),
              facing = "clockwise", niceFacing = T, cex = 0.25, adj = c(0, -0.2), font = 2)
}, track.height = 0.15, bg.border = NA)

##
pop.code <- cut(as.numeric(as.character(read.gdsn(index.gdsn(wheat, "samp.annot/Year")))),
                breaks = c(1800, 1900, 1920, 1940, 1960, 1980, 2000, 2020), labels = F)
pop.code[is.na(pop.code)] <- 8
pop.code <- factor(pop.code)
colours <- colors()[c(554, 53, 652, 48, 109, 26, 547, 153)]
# pie(rep(1,8), col = colours)
legend("topright", legend=c("Pre-1900", "1901-1920", "1921-1940", "1941-1960", "1961-1980", "1981-2000",
                            "2001-2016", "Unknown"), pch=19, col=colours, cex = 0.45, bg = colors()[468],
       title = "Period of Release")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:357.05, rep(0,358), 0.95:357.95, rep(1,358), border = colors()[468], lwd = 0.8,
              col = sapply(as.numeric(as.character(pop.code))[match(labels(neighbour), sample.id)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)

##
pop.code <- read.gdsn(index.gdsn(wheat, "samp.annot/BP"))
colours <- colors()[c(116, 554, 53, 26, 100, 12, 490, 42, 115, 139, 153, 451, 652, 114, 195, 557, 79, 48)]
legend("bottomright", legend=levels(pop.code), pch=19, col=colours, cex = 0.35, bg = colors()[142],
       title = "Breeding Program")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:357.05, rep(0,358), 0.95:357.95, rep(1,358), border = colors()[142], lwd = 0.8,
              col = sapply(as.numeric(pop.code)[match(labels(neighbour), sample.id)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)

##
pop.code <- read.gdsn(index.gdsn(wheat, "samp.annot/MC"))
colours <- colors()[c(195, 26, 49, 79, 652, 100, 419, 115, 114, 53, 554, 450, 153)]
legend("bottomleft", legend=levels(pop.code), pch=19, col=colours, cex = 0.45, bg = colors()[105],
       title = "Market Class")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:357.05, rep(0,358), 0.95:357.95, rep(1,358), border = colors()[105], lwd = 0.8,
              col = sapply(as.numeric(pop.code)[match(labels(neighbour), sample.id)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)

##
pop.code <- read.gdsn(index.gdsn(wheat, "samp.annot/designation"))
colours <- colours()[c(554, 115, 1, 26, 153, 450, 48, 114, 652, 153)]
legend("topleft", legend=levels(pop.code), pch=19, col=colours, cex = 0.45, bg = colors()[525],
       title = "Designation")
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:357.05, rep(0,358), 0.95:357.95, rep(1,358), border = colors()[525], lwd = 0.8,
              col = sapply(as.numeric(pop.code)[match(labels(neighbour), sample.id)], function(x) { return(colours[x]) } ))
}, track.height = 0.04, bg.border = NA)

##
load("Data\\Intermediate\\CCP\\full_kmeans.RDATA")
pop.code <- bests[match(labels(neighbour), sample.id)]
colours <-  colors()[c(53, 554, 100, 652, 114, 115)]
legend("top", title = "K-means Groups", legend=c("Group 1", "Group 2", "Group 3",
                                                 "Group 4", "Group 5", "Group 6"),
       pch=19, col=colours, cex = 0.4, bg = colors()[109], ncol = 2)
colours <-  colors()[c(554, 100, 114, 53, 652, 115)]
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(0.05:357.05, rep(0,358), 0.95:357.95, rep(1,358), border = colors()[109], lwd = 0.8,
              col = colours[pop.code])
}, track.height = 0.04, bg.border = NA)

max_height <- attr(neighbour, "height")
neighbour <- color_branches(neighbour, k = 5, col = c(brewer.pal(4, "Dark2"), "black"))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(neighbour, max_height = max_height)
}, track.height = 0.40, bg.border = NA)
legend("center", title = "NJ Groups", legend=c("Group A", "Group B",
                                               "Group C", "Group D"),
       pch=19, col=brewer.pal(4, "Dark2"), cex = 0.4)
# title("Fig 2: Neighbourjoining Dendrogram of\nWheat Varieties with Group Data", out = T, cex = 1.5, line = -2)
circos.clear()
dev.off()
snpgdsClose(wheat)


# library(ggtree)
# wheat <- snpgdsOpen("Data\\Formatted\\wheat_phys_subset_both.gds")
# sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))
# neighbour <- nj(1 - snpgdsIBS(wheat, autosome.only = F)$ibs)
# labels(neighbour) <- sample.id
# png("Analysis\\Figures\\dend\\full_nj.png", 
#     family="Times New Roman", width = 11, height = 11, pointsize = 20, units = "in", res = 500)
# ggtree(neighbour, layout="circular") + ggtitle("(Phylogram) circular layout")
# dev.off()
# snpgdsClose(wheat)