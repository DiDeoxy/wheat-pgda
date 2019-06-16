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
