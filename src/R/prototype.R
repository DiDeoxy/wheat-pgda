## ld decay

library(pgda)
library(SNPRelate)
library(tidyverse)
source(file.path("src", "R", "file_paths.R"))
# install.packages("sommer")
library(sommer)
library(mgcv)
library(pracma)

wheat_data <- snpgds_parse(phys_gds)

# Calcualte ld between all snps on each chromosome
wheat_gds <- snpgdsOpen(phys_gds)
ld_decay_data <- by(
  wheat_data$snp, wheat_data$snp$chrom, function (snp_data) {
    tibble(
      dist_mb = snp_data$pos %>% dist() %>% as.matrix() %>% as.vector() / 1e6,
      ld = snpgdsLDMat(
        wheat_gds, method = "composite", slide = -1, snp.id = snp_data$id
      )$LD %>% abs() %>% as.vector()
    )
  }
) %>% do.call(rbind, .)
snpgdsClose(wheat_gds)

model <- gam(ld ~ s(dist_mb, bs = "cs"), data = ld_decay_data)

predicted <- predict.gam(model, newdata = data.frame(dist_mb = 0:100))

ggplot() +
  geom_line(aes(seq_along(predicted), predicted))


# wheat_data <- snpgds_parse(phys_gds)

# # Calcualte ld between all snps on each chromosome
# wheat_gds <- snpgdsOpen(phys_gds)
# ld_decay_data <- by(
#   wheat_data$snp, wheat_data$snp$chrom, function (snp_data) {
#     neighbour_ld <- snpgdsLDMat(
#       wheat_gds, method = "composite", slide = -1, snp.id = snp_data$id
#     )$LD %>% Diag(., 1)
#     bad_markers <- -which(neighbour_ld > 0.90)
#     good_markers <- snp_data$id[bad_markers]
#     good_indices <- seq_along(snp_data$id)[bad_markers]
#     tibble(
#       dist_mb = snp_data$pos[good_indices] %>% dist() %>% as.matrix() %>%
#         as.vector() / 1e6,
#       ld = snpgdsLDMat(
#         wheat_gds, method = "composite", slide = -1, snp.id = good_markers
#       )$LD %>% abs() %>% as.vector()
#     )
#   }
# ) %>% do.call(rbind, .)
# snpgdsClose(wheat_gds)


# model <- gam(ld ~ s(dist_mb, bs = "cs"), data = ld_decay_data)

# predict.gam(model, newdata = data.frame(dist_mb = 0:100))


genome_length <- c(
    "1A" = 594102056, "1B" = 689851870, "1D" = 495451186,
    "2A" = 780798557, "2B" = 801256715, "2D" = 651852609,
    "3A" = 750843639, "3B" = 830829764, "3D" = 615552423,
    "4A" = 744588157, "4B" = 673617499, "4D" = 509857067,
    "5A" = 709773743, "5B" = 713149757, "5D" = 566080677,
    "6A" = 618079260, "6B" = 720988478, "6D" = 473592718,
    "7A" = 736706236, "7B" = 750620385, "7D" = 750620385
  ) %>% sum()


span_by_chrom(wheat_data$snp$chrom, wheat_data$snp$pos, diff = TRUE) %>% sum() / genome_length