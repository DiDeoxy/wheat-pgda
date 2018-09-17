library(tidyverse)
library(plyr)
library(GGally)
library(SNPRelate)
library(extrafont)

# load custom functions
source("src\\R_functions\\funcs_calc_stats.R")
source("src\\R_functions\\colour_sets.R")

# load the data from the gds object
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(data$pos, data$chrom, max)
max_genome_lengths <- c(max(chrom_lengths[seq(1, 19, 3)]), # A genome
                        max(chrom_lengths[seq(2, 20, 3)]), # B genome
                        max(chrom_lengths[seq(3, 21, 3)])) # D genome

# find the phi values of all markers in each of the three amova comparisons and
# add them to the data table
data_comps <- data
for (group in c("chrs_csws", "chrs_hrw", "csws_hrw")) {
  amova <- read_rds(str_c("Data\\Intermediate\\Adegenet\\", group,
                          "_amova.rds"))
  data_comps <- data_comps %>% add_column(!!(group) := phi_markers(amova))
}

# find the top 2.5% quantile for each of the three comparisons
signifs <- data_comps %>%
             summarise(chrs_csws = quantile(chrs_csws, prob = 0.975,
                                            na.rm = T),
                       chrs_hrw = quantile(chrs_hrw, prob = 0.975, na.rm = T),
                       csws_hrw = quantile(csws_hrw, prob = 0.975, na.rm = T))

# turn each value less than the top 2.5% quantile for the comparison to NA
data_comps$chrs_csws[which(data_comps$chrs_csws < signifs$chrs_csws)] <- NA
data_comps$chrs_hrw[which(data_comps$chrs_hrw < signifs$chrs_hrw)] <- NA
data_comps$csws_hrw[which(data_comps$csws_hrw < signifs$csws_hrw)] <- NA

# find which markers aren't na in each comparison
chrs_csws_markers <- which(! is.na(data_comps$chrs_csws))
chrs_hrw_markers <- which(! is.na(data_comps$chrs_hrw))
csws_hrw_markers <- which(! is.na(data_comps$csws_hrw))

# load the cluster data and make indexes of the individuals in each comparison
cluster <- read_rds("Data\\Intermediate\\dbscan\\full_hdbscan.rds")$cluster
hrs <- which(desig == "HRS")
index_chrs <- hrs[which(hrs %in% which(cluster == 5))]
sws <- which(desig == "SWS")
index_csws <- sws[which(sws %in% which(cluster == 3))]
index_hrw <- which(desig == "HRW")

# find the expected heterozygosity of each marker in each group and add to data
# we need one of each group in each of the comparisons, one is negated for 
# plotting below zero, like a mirror
data_eh <- data %>% 
             rbind(data) %>%
             add_column(eh_chrs_csws = c(eh(genotypes[, index_chrs]),
                                         -eh(genotypes[, index_csws]))) %>%
             add_column(eh_chrs_hrw = c(eh(genotypes[, index_chrs]),
                                        -eh(genotypes[, index_hrw]))) %>%
             add_column(eh_csws_hrw = c(eh(genotypes[, index_csws]), 
                                        -eh(genotypes[, index_hrw])))

# make all eh values not extreme in each comparison NA
data_eh$eh_chrs_csws[-c(chrs_csws_markers,
                        chrs_csws_markers + nrow(data))] <- NA
data_eh$eh_chrs_hrw[-c(chrs_hrw_markers,
                       chrs_hrw_markers + nrow(data))] <- NA
data_eh$eh_csws_hrw[-c(csws_hrw_markers,
                       csws_hrw_markers + nrow(data))] <- NA

# make the data long for easier plotting
data_long <- data_eh %>%
               gather(comparison, eh, c(eh_chrs_csws, eh_chrs_hrw, 
                                        eh_csws_hrw)) %>%
               arrange(chrom, pos)

# plot the EH values for each group in each comparison on each chromosome
lables <- c("CHRS vs CSWS", "CHRS vs HRW", "CSWS vs HRW")
legend_title <- "Comparisons"
plots <- by(data_long, data_long$chrom, function (data_chrom) {
  unique(data_chrom$comparison)
  chrom_num <- data_chrom$chrom[1]
  data_chrom %>%
    ggplot() +
      ylim(-0.5, 0.5) +
      xlim(0, max_genome_lengths[ifelse(chrom_num %% 3, chrom_num %% 3, 3)]) +
      geom_point(aes(pos, eh, colour = comparison, shape = comparison)) +
      geom_hline(yintercept = 0) +
      scale_colour_manual(legend_title,
                          labels = lables,
                          values = colours_comparisons_genes[1:3]) +
      scale_shape_manual(legend_title,
                         labels = lables,
                         values = points[1:3])
})

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "EH",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "EH Values of Markers in top 2.5% of AMOVA Phi Values by Comparison",
  legend = c(1, 1)
)

# plot the matrix
png(str_c("Results\\loci\\EH\\phi_eh.png"),
    family = "Times New Roman", width = 200, height = 287, pointsize = 5,
    units = "mm", res = 300)
plots_matrix + theme(legend.position = "bottom")
dev.off()