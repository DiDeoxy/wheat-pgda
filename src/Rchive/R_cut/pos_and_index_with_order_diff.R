source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
library(tidyverse)
library(GGally)
library(pgda)
library(circlize)

phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

snp_data <- tibble(
  phys_id = phys_data$snp$id, gen_id = gen_data$snp$id,
  chrom = phys_data$snp$chrom,
  phys_pos_mb = phys_data$snp$pos / 1e6,
  gen_pos_cm = gen_data$snp$pos / 100
)
snp_data_list <- split(snp_data, snp_data$chrom)

# calc the lengths of the different genomes and homoeologous sets
max_lengths_phys <- by(snp_data$phys_pos_mb, snp_data$chrom, max) %>% max_lengths()
max_lengths_gen <- by(snp_data$gen_pos_cm, snp_data$chrom, max) %>% max_lengths()
max_markers <- by(snp_data, snp_data$chrom, nrow) %>% max_lengths()

# create a function for making a gradient of colours
levels <- c(
  "0-4", "5-10", "11-15", "16-20", "21-50", "51-100", "101-500", "501-1500"
)
colour_levels <- colour_set[c(7, 4, 2, 3, 5, 1, 6, 19)]
names(colour_levels) <- levels

chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
  t() %>% as.vector()

circos.initialize(
  chroms, xlim = c(0, 1)
)
i <- 1
circos.track(
  ylim = c(0, 1), track.height = 0.15, bg.border = NA,
  panel.fun = function(x, y) {
    circos.text(
      1, 0, chroms[i], facing = "clockwise",
      niceFacing = TRUE, cex = 1, adj = c(0, -0.2), font = 2
    )
    i <<- i + 1
  }
)
i <- 1
circos.track(
  track.height = 0.3,
  ylim = c(
    0,
    by(snp_data$gen_pos_cm, snp_data$chrom, max) %>%
      as.list() %>% unlist() %>% max()
  ),
  panel.fun = function (x, y) {
    gen_to_phys_order <- match(
      snp_data_list[[i]]$phys_id, snp_data_list[[i]]$gen_id
    )
    order_diff <- (gen_to_phys_order - 1:length(gen_to_phys_order)) %>% abs()

    cols <- colour_levels[
      (
        cut(order_diff, c(-1, 4, 10, 15, 20, 50, 100, 500, 1500), levels ) %>%
          as.integer()
      )
    ]
    circos.points(
      x = snp_data_list[[i]]$phys_pos_mb / max(snp_data_list[[i]]$phys_pos_mb),
      y = snp_data_list[[i]]$gen_pos_cm[gen_to_phys_order], col = cols,
      pch = 19, cex = 0.5
    )
    i <<- i + 1
  }
)
i <- 1
circos.track(
  track.height = 0.3,
  ylim = c(
    0, by(snp_data, snp_data$chrom, nrow) %>% as.list() %>% unlist() %>% max()
  ),
  panel.fun = function (x, y) {
    gen_to_phys_order <- match(
      snp_data_list[[i]]$phys_id, snp_data_list[[i]]$gen_id
    )
    order_diff <- (gen_to_phys_order - 1:length(gen_to_phys_order)) %>% abs()

    cols <- colour_levels[
      (
        cut(order_diff, c(-1, 4, 10, 15, 20, 50, 100, 500, 1500), levels ) %>%
          as.integer()
      )
    ]
    circos.points(
      x = 1:length(gen_to_phys_order) / length(gen_to_phys_order),
      y = gen_to_phys_order, col = cols, pch = 19, cex = 0.5
    )
    i <<- i + 1
  }
)
circos.clear()

