library(plyr)
library(tidyverse)
library(GGally)
library(ggrepel)
library(extrafont)
library(SNPRelate)

# load custom functions
source("src\\R_functions\\funcs_calc_stats.R")
source("src\\R_functions\\colour_sets.R")

# size of a megabase, used to divide the bp positions
mb <- 1000000

# load the data from the gds object
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")
snp_pos <- snp_pos / mb

# create a data frame of all the 
data <- tibble(id = snp_id, chrom = snp_chrom, pos = snp_pos)

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(data$pos, data$chrom, max)
max_genome_lengths <- data.frame(A = max(chrom_lengths[seq(1, 19, 3)]),
                                 B = max(chrom_lengths[seq(2, 20, 3)]),
                                 D = max(chrom_lengths[seq(3, 21, 3)]))

# find the phi values of all markers in each of the three amova comparisons and
# add them to the data table
for (group in c("chrs_csws", "chrs_hrw", "csws_hrw")) {
  amova <- read_rds(str_c("Data\\Intermediate\\Adegenet\\", group,
                          "_amova.rds"))
  data <- data %>% add_column(!!(group) := phi_markers(amova))
}

# find the top 2.5% quantile for each of the three comparisons
signifs <- data %>%
             summarise(chrs_csws = quantile(chrs_csws, prob = 0.975,
                                            na.rm = T),
                       chrs_hrw = quantile(chrs_hrw, prob = 0.975, na.rm = T),
                       csws_hrw = quantile(csws_hrw, prob = 0.975, na.rm = T))

# turn each value less than the top 2.5% quantile for the comparison to NA
data$chrs_csws[which(data$chrs_csws < signifs$chrs_csws)] <- NA
data$chrs_hrw[which(data$chrs_hrw < signifs$chrs_hrw)] <- NA
data$csws_hrw[which(data$csws_hrw < signifs$csws_hrw)] <- NA

# load gene data
main_genes <- read_rds(
  "Data\\Intermediate\\Aligned_genes\\top_main_genes_selected_corrected.rds")
resi_genes <- read_rds(
  str_c("Data\\Intermediate\\Aligned_genes\\top_resistance_genes_selected",
        "_corrected.rds")
)

# add min phi values of plotting to each gene so they appear at bottom of plots
main_genes <- cbind(main_genes, main_genes = min(signifs))
resi_genes <- cbind(resi_genes, resi_genes = min(signifs))

# make the data set long for easier plotting
data_long <- data %>%
             gather(comparison, phi, c(chrs_csws, chrs_hrw, csws_hrw))

# add the genes in and make longer
data_long_genes <- data_long %>%
                   rbind.fill(main_genes) %>%
                   rbind.fill(resi_genes) %>%
                   arrange(chrom, comparison, pos) %>%
                   gather(gene_type, phi2, c(main_genes, resi_genes)) %>%
                   as.tibble()

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
chrom_num <- 0
lables <- c("CHRS vs CSWS", "CHRS vs HRW", "CSWS vs HRW", "Phenotype Genes",
            "Resistance Genes")
legend_title <- "Comparisons and Genes"
plots <- by(data_long_genes, data_long_genes$chrom, function (data_chrom) {
  chrom_num <<- data_chrom$chrom[1]
  data_chrom %>%
    ggplot() +
      ylim(min(signifs), 1) +
      xlim(0, max_genome_lengths[ifelse(chrom_num %% 3, chrom_num %% 3, 3)]) +
      geom_point(aes(pos, phi, colour = comparison, shape = comparison),
                 size = 0.5) +
      geom_point(aes(pos, phi2, colour = gene_type, shape = gene_type)) +
      geom_text_repel(aes(pos, phi2, colour = gene_type, label = id),
                      angle = 90, hjust = 0, vjust = 1, size = 3,
                      segment.colour = "black", nudge_y = 0.07,
                      show.legend = FALSE) +
      scale_colour_manual(legend_title,
                          labels = lables,
                          values = colours_comparisons_genes) +
      scale_shape_manual(legend_title,
                        labels = lables,
                        values = c(15, 16, 18, 17, 8))
})

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "Phi",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Markers in Top 2.5% of AMOVA Phi Values By Comparison",
  legend = c(1, 1)
)

# plot the matrix
png(str_c("Results\\loci\\amova\\full_amova.png"),
    family = "Times New Roman", width = 200, height = 287, pointsize = 5,
    units = "mm", res = 300)
plots_matrix + theme(legend.position = "bottom")
dev.off()

# # find closest markers to genes
# data_omit_genes <- data_long %>%
#                    na.omit() %>%
#                    rbind.fill(main_genes) %>%
#                    rbind.fill(resi_genes) %>%
#                    select(id, comparison, chrom, pos) %>%
#                    arrange(chrom, comparison, pos) %>%
#                    as.tibble()

# by(data_omit_genes, data_omit_genes$chrom, function (markers_genes) {

#   if ("YR10" %in% markers_genes$id) {
#     closest = which.min(abs(
#                 markers_genes$pos[which(! (is.na(marker_genes$comparison))] - 
#                 markers_genes$pos[which(markers_genes$id == "YR10")]))
#     return(markers_genes$id[closest])
#   }
# })

# summary of phi values
# phi_table <- function () {
#   snps[[group]] <- by(data.frame(snp_id, phi_all, snp_pos), snp_chrom, function (chr) {
#     return(chr[which(chr[,2] > signif),3])
#   })

#   phiStats <- by(phi_all, snp_chrom, function (chr) {
#     return(cbind(mean(abs(diff(chr))), sd(abs(diff(chr)))))
#   })
#   order <- c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3), 22)
#   labels <- c(as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep=""))), "All")
#   phiTable <- data.frame()
#   phiTable <- do.call("rbind", phiStats)
#   phiTable <- rbind(phiTable, cbind(mean(phi_all), sd(phi_all)))
#   row.names(phiTable) <- labels
#   colnames(phiTable) <- c("Mean", "Standard Deviation")
#   phiTable <- phiTable[order,]
#   write.csv(phiTable, file = paste0("Data\\Intermediate\\phi stats\\", group, ".csv"))
# }
