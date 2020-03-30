# load the needed functions and file paths
source(file.path("src", "R", "file_paths.R"))
import::from(
  circlize, "circos.initialize", "circos.clear", "circos.dendrogram", 
  "circos.par", "circos.rect", "circos.text", "circos.track"
)
import::from(dendextend, "color_branches", "set")
import::from(magrittr, "%>%")
import::from(stringi, "stri_rand_strings")
import::from(stringr, "str_c")

################################################################################
# set parameters

path <- getwd() # change to wherever you want to the figure

# numbers of things
n_indivs <- sample(100:500, 1)
n_markers <- sample(seq(1000, 5000, 100), 1)
class_size <- sample(3:12, sample(1:4, 1))
n_dend_clusts <- sample(1:15, 1)

# colours
colour_set <- c(
    "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
    "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080",
    "#FFFFFF", "#000000"
)
# pie(rep(1, 22), col = colour_set)

class_colours <- list()
for (i in seq_along(class_size)) {
  class_colours[[i]] <- sample(colour_set, class_size[i])
}
dend_colours <- sample(colour_set, n_dend_clusts)
bg_colours <- sample(colour_set, length(class_size))

legend_position <- c(
  "topright", "bottomright", "bottomleft", "bottomright"
)

################################################################################
# make a mock data set

# really bad genotype simulation
genotypes <- sample(c(0, 2), n_markers * n_indivs, replace = TRUE) %>%
  matrix(nrow = n_markers, ncol = n_indivs)

indivs <- stri_rand_strings(n_indivs, 10, pattern = "[A-Za-z]")
dend_clusts <- sample(
  stri_rand_strings(n_dend_clusts, 5, pattern = "[A-Za-z]"), n_indivs,
  replace = TRUE
) %>% as.factor()

# simulating different classifications of indivs
class_data <- list()
for (i in seq_along(class_size)) {
  class_data[[i]] <- sample(
    stri_rand_strings(class_size[i], 5, pattern = "[A-Za-z]"), n_indivs,
    replace = TRUE
  ) %>% as.factor()
}

################################################################################

############### comment out everything above to freeze data set ################

################################################################################
# make the dendrogram

# calcualte the fraction of identity by state (IBS) of markers between each
# individual
ibs_mat <- matrix(NA, nrow = n_indivs, ncol = n_indivs)
for (i in 1:ncol(genotypes)) {
  if (i < ncol(genotypes)) {
    for (j in (i + 1):ncol(genotypes)) {
      ibs_mat[j, i] <- sum(genotypes[, i] == genotypes[, j]) / nrow(genotypes)
    }
  }
}
ibs_dist <- as.dist(ibs_mat)

# create the dendrogram from the IBS matrix
upgma_dend <- hclust(ibs_dist) %>%
  as.dendrogram(method = "average") %>%
  color_branches(k = n_dend_clusts, col = dend_colours) %>%
  set("branches_lwd", 1.5)

# need this order so that the correct labels can be applied on the plot
label_order <- order.dendrogram(upgma_dend)

################################################################################
# function for drawing a track of rectangles

draw_rects <- function(pop_code, colour_subset, label_order, border) {
  leng <- length(label_order)
  # defines how the recatangles should be sized and what colours they should be
  rects <- function(x, y) {
    circos.rect(
      0.05:(leng - 0.95), rep(0, leng),
      0.95:(leng - 0.05), rep(1, leng),
      border = border,
      lwd = 0.9,
      col = sapply(
        as.numeric(pop_code)[label_order],
        function(z) {
          return(colour_subset[z])
        }
      )
    )
  }

  # actually draws the track
  circos.track(
    ylim = c(0, 1),
    panel.fun = rects(x, y),
    track.height = 0.04,
    bg.border = NA
  )
}

################################################################################
# draw the circos plot

png(
  file.path(path, "dend.png"), family = "Times New Roman", width = 210,
  height = 210, pointsize = 15, units = "mm", res = 192
)
circos.par(
  cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, track.margin = c(0.005, 0.005)
)
circos.initialize("foo", xlim = c(0, length(label_order)), sector.width = 1)

# plot the names of each individual at the correct position
circos.track(
  ylim = c(0, 1), track.height = 0.15, bg.border = NA,
  panel.fun = function(x, y) {
    circos.text(
      1:length(label_order), rep(0, length(label_order)),
      indivs[label_order], facing = "clockwise",
      niceFacing = TRUE, cex = 0.25, adj = c(0, -0.2), font = 2
    )
  }
)

# draw each of the sets of class data as a separate track of rectangles
for (i in seq_along(class_size)) {
  draw_rects(
    class_data[[i]], class_colours[[i]], label_order, bg_colours[i]
  )
}

# draw the dendrogram
max_height <- max(attr(upgma_dend, "height"))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(upgma_dend, max_height = max_height)
}, track.height = 0.3, bg.border = NA)

circos.clear()

# draw a legend for each class
bg_colours <- sample(colour_set, length(class_size))
for (i in seq_along(class_size)) {
  legend(
    legend_position[i], legend = levels(class_data[[i]]), box.lwd = 2,
    title = str_c("Class ", i), pch = 19, col = class_colours[[i]], cex = 0.45,
    bg = bg_colours[i]
  )
}

legend(
  "center", legend = levels(dend_clusts), box.lwd = 2,
  title = str_c("UPGMA Clusters"), pch = 19, col = dend_colours,
  cex = 0.45
)

dev.off()
