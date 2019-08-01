colour_set <- c(
    "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
    "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080",
    "#FFFFFF", "#000000"
)
# png(file.path("results", "colours.png"))
# pie(rep(1, 22), col = colour_set)
# dev.off()

# defining the colour sets
# colours_order_diff <- colour_set[c(7, 4, 2, 5, 1, 6, 19, 22)]
colours_intervals <- colour_set[c(2, 4, 5, 1, 6, 19, 13, 15)]
colours_era <- colour_set[c(1, 5, 3, 2, 4, 6, 22)]
colours_bp <- colour_set[c(3, 1, 5, 4, 12, 2, 20, 22, 6, 7, 8, 9, 13)]
colours_mc <- colour_set[c(5, 4, 3, 6, 2, 11, 1, 7, 22)]
colours_pheno <- colour_set[c(1, 2, 7, 16, 12, 4, 9, 22)]
colours_pheno_ar <- colour_set[c(1, 2, 7, 12, 4, 9, 22)]
# colours_dend <- colour_set[c(22, 4, 22, 1, 3, 5, 1, 2, 22)]
colours_dend <- colour_set[c(5, 1, 4, 22, 22, 2, 1, 22, 3)]
colours_hdbscan <- colour_set[c(22, 2, 4, 3, 5, 1)]
# clusters 2 and 3 are swapped on the legend of the dendrogram for writeup
# purposes
colours_hdbscan_legend <- colour_set[c(22, 2, 4, 3, 5, 1)]
colours_hdbscan_ar_legend <- colour_set[c(2, 4, 3, 5, 1, 22)]
colours_chroms <- colour_set[c(1, 5, 3, 2, 4, 6, 8)]
colours_groups_genes <- colour_set[c(1, 2, 4, 15, 19)]
colours_comps_genes <- colour_set[c(1, 2, 4, 20, 15, 19)]
# colours_comps_genes <- colour_set[c(5, 1, 17, 2, 11, 4, 20, 15, 19)]
# pie(rep(1, 9), col = colours_dend)