source("wheat-pgda/src/R/file_paths.R")
library(ggplot2)
import::from(ggrepel, "geom_label_repel")
import::from(gridExtra, "grid.arrange")
import::from(magrittr, "%>%")
import::from(parallel, "mclapply")
import::from(pgda, "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(readr, "read_csv")
import::from(SNPRelate, "snpgdsLDMat", "snpgdsClose", "snpgdsOpen")
import::from(gdsfmt, "index.gdsn", "read.gdsn")

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

calc_pseudo_gen_pos <- function (ld_mat, cores) {
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

sparsity <- lapply(ld_mats, calc_pseudo_gen_pos, 8)

gen_pos <- lapply(sparsity, cumsum)

phys_data <- snpgds_parse(phys_gds)

max_phys_lengths <- span_by_chrom(
  phys_data$snp$chrom, phys_data$snp$pos
) %>% max_lengths() / 1e6

landmarks <- read_csv(file.path(intermediate, "centromeres.csv")) %>% cbind(id = "Centromere", base = 0.5) %>% split(.$chrom)

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
  family = "Times New Roman", width = 900, height = 900, pointsize = 5,
  units = "mm", res = 192
)
grid.arrange(
  grobs = plots, nrow = 7, ncol = 3
)
dev.off()