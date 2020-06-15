library(SNPRelate)
library(scrime)
library(ape)
library(ggtree)
library(extrafont)
library(dendextend)

## loading the population data
wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset_imputed.gds")
genotypes <- read.gdsn(index.gdsn(wheat.subset, "genotype"))
sample.id <- as.character(read.gdsn(index.gdsn(wheat.subset, "sample.id")))
snp.pos <- read.gdsn(index.gdsn(wheat.subset, "snp.position"))
snp.chr <- read.gdsn(index.gdsn(wheat.subset, "snp.chromosome"))
meta.data <- cbind(BP = as.character(read.gdsn(index.gdsn(wheat.subset, "samp.annot/BP"))),
                   Year = as.character(read.gdsn(index.gdsn(wheat.subset, "samp.annot/Year"))),
                   MC = as.character(read.gdsn(index.gdsn(wheat.subset, "samp.annot/MC"))),
                   Desig = as.character(read.gdsn(index.gdsn(wheat.subset, "samp.annot/designation"))))
rownames(meta.data) <- sample.id

## Clustering
# make the dist object
wheat.dist <- as.dist(1-snpgdsIBS(wheat.subset)$ibs)
snpgdsClose(wheat.subset)

# reformat some stuff
desig <- split(sample.id, meta.data[,3])

# neighbour joining
nj <- nj(wheat.dist)
labels(nj) <- sample.id
nj <- groupOTU(nj, desig)

## Colour it
colours <- colors()[c(554, 115, 49, 26, 153, 450, 114, 145, 153)]
colours <- colors()[c(195, 26, 49, 79, 145, 100, 419, 115, 114, 53, 554, 450, 153)]
names(colours) <- names(desig)

## plot it
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\dend\\NJ_desig.png", 
    family="Times New Roman", width = 11, height = 11, pointsize = 20, units = "in", res = 300)
ggtree(nj, aes(colour=group), layout = "circular", branch.length = "none") + 
  geom_tiplab2(size = 2) + scale_colour_manual(name = "Designations", values=colours) +
  theme(legend.position="right")
# gheatmap(bhc.phylo, meta.data, offset = 2, width=0.5)
dev.off()

# genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
# genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
# class(meta.data)
# class(genotype)
