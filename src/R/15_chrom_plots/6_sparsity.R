source("wheat-pgda/src/R/file_paths.R")
source("wheat-pgda/src/R/colour_sets.R")
library(ggplot2)
import::from(gdsfmt, "index.gdsn", "read.gdsn")
import::from(ggalt, "geom_xspline")
import::from(ggrepel, "geom_label_repel")
import::from(gridExtra, "grid.arrange")
import::from(magrittr, "%>%")
import::from(parallel, "mclapply")
import::from(pgda, "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(readr, "read_csv")
import::from(SNPRelate, "snpgdsLDMat", "snpgdsClose", "snpgdsOpen")
import::from(tibble, "tibble")

phys_data <- snpgdsOpen(phys_gds)

chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>% 
  t() %>% as.vector()

chrom <- as.character(read.gdsn(index.gdsn(phys_data, "snp.chromosome")))

for (i in 1:21) {
  chrom[which(chrom == i)] <- chroms[i]
}

var_sets <- as.character(read.gdsn(index.gdsn(phys_data, "snp.id"))) %>%
  split(chrom)
phys_pos <- as.numeric(read.gdsn(index.gdsn(phys_data, "snp.position"))) %>%
  split(chrom)

pane <- 10
window <- pane * 2 + 1

ld_mats <- lapply(var_sets, function (var_set) {
  snpgdsLDMat(phys_data, slide = window, snp.id = var_set, num.thread = 8)$LD
})

snpgdsClose(phys_data)

calc_sparsity <- function (ld_mat, cores) {
  int_window <- min(window, ncol(ld_mat))
  mclapply(1:ncol(ld_mat), function (i) {
    if (i - pane >= 0 && i + pane <= ncol(ld_mat)) {
      cols <- c((i - pane):i, rep(i, pane - 1))
      rows <- c(pane:1, 1:(pane))
    } else if (i - pane < 0) {
      cols <- 1:(int_window - 1)
      if (i == 1) {
        rows <- c(1:length(cols))
      } else {
        rows <- c(sum(cols < i):1, 1:sum(cols >= i))
      }
      cols[which(cols > i)] <- i
    } else {
      cols <- (
          ncol(ld_mat) - (int_window - 1)
      ):(
          ncol(ld_mat) - 1
      )
      if (i == ncol(ld_mat)) {
        rows <- c(length(cols):1)
      } else {
        rows <- c(sum(cols < i):1, 1:sum(cols >= i))
      }
      cols[which(cols > i)] <- i
    }
    indices <- cbind(rows, cols)
    lds <- abs(ld_mat[indices])
    lds[is.nan(lds)] <- NA
    sparsity <- -log(mean(lds, na.rm = TRUE))
    if (is.nan(sparsity)) {
        return(0)
    } else {
        return(sparsity)
    }
  }, mc.cores = cores) %>% unlist() %>% as.numeric()
}

sparsity <- lapply(ld_mats, calc_sparsity, 8)

phys_data <- snpgds_parse(phys_gds)

max_phys_lengths <- span_by_chrom(
  phys_data$snp$chrom, phys_data$snp$pos
) %>% max_lengths() / 1e6

landmarks <- read_csv(file.path(intermediate, "centromeres.csv")) %>%
  cbind(id = "Centromere", base = 0.5) %>%
  split(.$chrom)

sparsity_by_chrom <- tibble(chrom = phys_data$snp$chrom, pos_mb = phys_data$snp$pos / 1e6, sparsity = sparsity %>% unlist())

################################################################################

plots <- lapply(chroms, function (chrom) {

  max_pos_mb <- max_phys_lengths[[
    ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
  ]]

  ggplot() +
    geom_point(aes(phys_pos[[chrom]] / 1e6, sparsity[[chrom]])) +
    labs(y = "Genetic Sparsity", x = "Pseudo-Chromosome Position", title = chrom) +
    expand_limits(x = c(0, max_pos_mb)) +
    scale_x_continuous(
      breaks = seq(0, max_pos_mb, by = 25), expand = c(0.01, 0.01)
    ) +
    geom_vline(aes(xintercept = landmarks[[chrom]]$pos_mb)) +
    geom_label_repel(
      aes(
        landmarks[[chrom]]$pos_mb,
        landmarks[[chrom]]$base * max(sparsity[[chrom]]),
        label = landmarks[[chrom]]$id
      )
    )
})

png(
   "/workspace/results/sparsity.png",
  family = "Times New Roman", width = 2480, height = 3508
)
grid.arrange(
  grobs = plots, nrow = 7, ncol = 3
)
dev.off()

# ################################################################################

# groups <- lapply(phys_data$snp$chrom, function (chrom) {
#   ifelse(grepl("1", chrom), "1",
#       ifelse(grepl("2", chrom), "2",
#         ifelse(grepl("3", chrom), "3",
#           ifelse(grepl("4", chrom), "4",
#             ifelse(grepl("5", chrom), "5",
#               ifelse(grepl("6", chrom), "6", "7")
#             )
#           )
#         )
#       )
#     )
# }) %>% unlist()

# plots <- by(sparsity_by_chrom, groups, function (group) {
#   group2 <- by(group, group$chrom, function (chrom) {
#     chrom_spline <- smooth.spline(chrom$pos_mb, chrom$sparsity, spar = 0.9)
#     tibble(chrom = chrom$chrom[1], x = chrom_spline$x, y = chrom_spline$y)
#   }) %>% as.list() %>% do.call(rbind, .)

#   ggplot() +
#     geom_point(data = group, aes(pos_mb, sparsity, colour = chrom), alpha = 0.1) +
#     geom_line(data = group2, aes(x, y, color = chrom), size = 2) +
#     scale_colour_manual(
#       values = colour_set[c(1, 2, 4)]
#     )
# }) %>% as.list()

# png(
#    "/workspace/results/sparsity.png",
#   family = "Times New Roman", width = 1240, height = 1754
# )
# grid.arrange(
#   grobs = plots, nrow = 7, ncol = 1
# )
# dev.off()

# ################################################################################

# genomes <- lapply(phys_data$snp$chrom, function (chrom) {
#   ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
# }) %>% unlist()

# plots <- by(sparsity_by_chrom, genomes, function (genome) {
#   genome2 <- by(genome, genome$chrom, function (chrom) {
#     chrom_spline <- smooth.spline(chrom$pos_mb, chrom$sparsity, spar = 0.5)
#     tibble(chrom = chrom$chrom[1], x = chrom_spline$x, y = chrom_spline$y)
#   }) %>% as.list() %>% do.call(rbind, .)

#   ggplot() +
#     geom_point(data = genome, aes(pos_mb, sparsity, colour = chrom), alpha = 0.1) +
#     geom_line(data = genome2, aes(x, y, color = chrom), size = 1)
# }) %>% as.list()

# png(
#    "/workspace/results/sparsity.png",
#   family = "Times New Roman", width = 1240, height = 1754
# )
# grid.arrange(
#   grobs = plots, nrow = 3, ncol = 1
# )
# dev.off()