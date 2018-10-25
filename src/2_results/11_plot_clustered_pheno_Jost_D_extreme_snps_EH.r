library(tidyverse)
library(plyr)
library(GGally)
library(SNPRelate)
library(extrafont)

# load custom functions
source("src\\R_functions\\funcs_gds_parse_create.R")
source("src\\R_functions\\funcs_calc_stats.R")
source("src\\R_functions\\colour_sets.R")

# load the data from the gds object
wheat_data <- parse_gds("phys_subset_sample")

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(wheat_data$snp$pos, wheat_data$snp$chrom, max)
max_genome_lengths <- c(max(chrom_lengths[seq(1, 19, 3)]), # A genome
                        max(chrom_lengths[seq(2, 20, 3)]), # B genome
                        max(chrom_lengths[seq(3, 21, 3)])) # D genome

# add columns of the D values of each comparison
for (group in c("chrs_csws", "chrs_chrw", "csws_chrw")) {
  comp_Jost_D <- read_rds(str_c("Data\\Intermediate\\mmod\\", group,
    "_Jost_D.rds"))[[1]]
  wheat_data$snp <- wheat_data$snp %>%
    add_column(!!str_c(group, "_D") := comp_Jost_D)
}

# find the top 2.5% quantile for D for each of the three comparisons
extremes <- wheat_data$snp %>%
  summarise(
    chrs_csws_D = quantile(chrs_csws_D, prob = 0.975, na.rm = T),
    chrs_chrw_D = quantile(chrs_chrw_D, prob = 0.975, na.rm = T),
    csws_chrw_D = quantile(csws_chrw_D, prob = 0.975, na.rm = T)
  )

# load the cluster data and make indexes of the individuals in each comparison
cluster <- read_rds("Data\\Intermediate\\hdbscan\\wheat_hdbscan.rds")$cluster
pheno_indices <- list()
pheno_indices[["CHRW"]] <- which(
  wheat_data$sample$annot$pheno == "HRW" & cluster == 1
)
pheno_indices[["CSWS"]] <- which(
  wheat_data$sample$annot$pheno == "SWS" & cluster == 2
)
pheno_indices[["CHRS"]] <- which(
  wheat_data$sample$annot$pheno == "HRS" & cluster == 5
)

# find the expected heterozygosity of each marker in each group and add to data
# we need one of each group in each of the comparisons, one is negated for 
# plotting below zero, like a mirror
wheat_data$snp <- rbind(wheat_data$snp, wheat_data$snp) %>%
  add_column(
    chrs_csws_eh = c(
      calc_eh(wheat_data$genotypes[, pheno_indices$CHRS]),
      -calc_eh(wheat_data$genotypes[, pheno_indices$CSWS])
    )
  ) %>%
  add_column(
    chrs_chrw_eh = c(
      calc_eh(wheat_data$genotypes[, pheno_indices$CHRS]),
      -calc_eh(wheat_data$genotypes[, pheno_indices$CHRW])
    )
  ) %>%
  add_column(
    csws_chrw_eh = c(
      calc_eh(wheat_data$genotypes[, pheno_indices$CSWS]),
      -calc_eh(wheat_data$genotypes[, pheno_indices$CHRW])
    )
  ) %>%
  as.tibble()

# make all eh values that dont have extreme D values in each comparison NA
wheat_data$snp$chrs_csws_eh[
  which(wheat_data$snp$chrs_csws_D < extremes$chrs_csws_D)] <- NA
wheat_data$snp$chrs_chrw_eh[
  which(wheat_data$snp$chrs_chrw_D < extremes$chrs_chrw_D)] <- NA
wheat_data$snp$csws_chrw_eh[
  which(wheat_data$snp$csws_chrw_D < extremes$csws_chrw_D)] <- NA

mins <- c(
  chrs_csws_eh = min(wheat_data$snp$chrs_csws_eh),
  chrs_chrw_eh = min(wheat_data$snp$chrs_chrw_eh),
  csws_chrw_eh = min(wheat_data$snp$csws_chrw_eh)
)

# make the data long for easier plotting
wheat_data$snp_long <- wheat_data$snp %>%
  gather(
    comparison, eh, c(chrs_csws_eh, chrs_chrw_eh, csws_chrw_eh)
  ) %>%
  arrange(chrom, pos)

# plot the EH values for each group in each comparison on each chromosome
lables <- c("CHRS vs CHRW", "CHRS vs CSWS", "CSWS vs CHRW")
legend_title <- "Comparisons"
plots <- by(wheat_data$snp_long, wheat_data$snp_long$chrom,
  function (data_chrom) {
    chrom_num <- data_chrom$chrom[1]
    data_chrom %>%
      ggplot() +
        ylim(-0.5, 0.5) +
        xlim(
          0,
          max_genome_lengths[ifelse(chrom_num %% 3, chrom_num %% 3, 3)]
        ) +
        geom_point(aes(pos, eh, colour = comparison, shape = comparison),
          size = 0.75
        ) +
        geom_hline(yintercept = 0) +
        scale_colour_manual(
          legend_title, labels = lables,
          values = colours_comparisons_genes[1:3]
        ) +
        scale_shape_manual(
          legend_title, labels = lables, values = points[1:3]
        )
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "EH",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "EH Values of Markers in Top 2.5% of Jost D Values by Comparison",
  legend = c(1, 1)
)

# plot the matrix
png(str_c("Results\\loci\\EH\\extreme_D_eh_comps.png"),
    family = "Times New Roman", width = 210, height = 267, pointsize = 5,
    units = "mm", res = 300)
plots_matrix + theme(legend.position = "bottom")
dev.off()