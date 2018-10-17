library(SNPRelate)
library(ggplot2)
library(GGally)

# function for creating ld scatter plots
create_plots <- function(full, marker_data) {
  # chr1 <- which(marker_data$chrom == 1)
  by(marker_data, marker_data$chrom, function(chrom) {
    ld_mat <- abs(snpgdsLDMat(
      full,
      method = "composite", snp.id = chrom$id, slide = -1
    )$LD)

    dist_ld <- unlist(sapply(2:(nrow(chrom)), function(i) {
      sapply(1:(i - 1), function(j) {
        c(chrom$pos[i] - chrom$pos[j], ld_mat[i, j])
      })
    }))
    dist_ld <- tibble(
      distance = dist_ld[seq(1, length(dist_ld), 2)],
      ld = dist_ld[seq(2, length(dist_ld), 2)]
    )

    set.seed(1000)
    rows <- sample(nrow(dist_ld), nrow(dist_ld))
    ggplot(dist_ld[rows, ], aes(distance, ld)) +
      geom_point(size = 0.10) +
      geom_density_2d()
  })
}

gds <- "Data\\Intermediate\\GDS\\full_phys.gds"
source("src\\R_functions\\data_loading.R")
full <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys.gds")
# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  create_plots(full, marker_data), nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = "LD", xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Pairwise LD between all Markers on each Chrom by Distance"
  # legend = c(1, 1)
)
# plot the matrix
png(str_c("Results\\loci\\LD\\full_LD_decay.png"),
    family = "Times New Roman", width = 200, height = 287, pointsize = 5,
    units = "mm", res = 300)
plots_matrix
dev.off()
snpgdsClose(full)


gds <- "Data\\Intermediate\\GDS\\full_phys_floor_pruned.gds"
source("src\\R_functions\\data_loading.R")
full <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys_floor_pruned.gds")
# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  create_plots(full, marker_data), nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = "LD", xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Pairwise LD between all Markers on each Chrom by Distance"
  # legend = c(1, 1)
)
# plot the matrix
png(str_c("Results\\loci\\LD\\floor_LD_decay.png"),
    family = "Times New Roman", width = 200, height = 287, pointsize = 5,
    units = "mm", res = 300)
plots_matrix
dev.off()
snpgdsClose(full)