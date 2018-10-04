colour_set <- c(
    "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
    "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080",
    "#FFFFFF", "#000000"
)
# pie(rep(1, 22), col = colour_set)

# defining the colour sets
colours_era <- colour_set[c(1, 5, 3, 2, 4, 6, 22)]
colours_bp <- colour_set[c(3, 1, 5, 7, 12, 2, 20, 22, 6, 4, 8, 9, 13)]
colours_mc <- colour_set[c(5, 4, 11, 13, 2, 6, 15, 8, 9, 10, 1, 12, 22)]
colours_pheno <- colour_set[c(1, 6, 4, 16, 12, 7, 9, 22)]
colours_dend <- colour_set[c(11, 1, 7, 15, 6, 5, 17)]
colours_hdbscan <- colour_set[c(22, 6, 7, 5, 17, 1)]
# clusters 2 and 3 are swapped on the legend of the dendrogram for writeup
# purposes
colours_hdbscan_legend <- colour_set[c(22, 6, 7, 5, 17, 1)]
colours_chroms <- colour_set[c(1, 5, 3, 2, 4, 6, 8)]
colours_comparisons_genes <- colour_set[c(1, 2, 4, 15, 19)]
points <- c(15, 16, 18, 17, 8)
pie(rep(1, 7), col = colours_dend)