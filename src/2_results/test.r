library(tidyverse)
library(plyr)
library(GGally)
library(ggrepel)
library(extrafont)
library(SNPRelate)

# load custom functions
source("src/R_functions/funcs_calc_stats.R")
source("src/R_functions/funcs_gds_parse_create.R")

# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(wheat_data$snp$pos_mb, wheat_data$snp$chrom, max)
max_genome_lengths <- c(
  A = max(chrom_lengths[seq(1, 19, 3)]),
  B = max(chrom_lengths[seq(2, 20, 3)]),
  D = max(chrom_lengths[seq(3, 21, 3)])
)

# find the phi values of each marker in each Gene and add to data set
for (group in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
  gene_Jost_D <- read_rds(str_c("Data/Intermediate/mmod/", group,
    "_Jost_D.rds"))[[1]]
  wheat_data$snp <- wheat_data$snp %>% add_column(!!group := gene_Jost_D)
}

# find the significance threshold for each Gene
signifs <- wheat_data$snp %>%
  summarise(
    Lr34 = quantile(Lr34, prob = 0.975, na.rm = TRUE),
    Lr22a = quantile(Lr22a, prob = 0.975, na.rm = TRUE),
    Lr21 = quantile(Lr21, prob = 0.975, na.rm = TRUE),
    Lr10 = quantile(Lr10, prob = 0.975, na.rm = TRUE),
    Lr1 = quantile(Lr1, prob = 0.975, na.rm = TRUE)
  )

# gene_signif_freq <- by(wheat_data$snp, wheat_data$snp$chrom,
#   function(snp_data) {
#     ret <- tibble(
#       gene = character(), chrom = integer(), pos = double(), freq = double()
#     )
#     for (gene in names(signifs)) {
#       for (i in 5:(length(snp_data$id) - 5)) {
#         ret <- ret %>%
#           add_row(
#             gene = gene, chrom = snp_data$chrom[1], 
#             pos = mean(snp_data$pos_mb[(i - 4):(i + 5)]),
#             freq = sum(snp_data[[gene]][(i - 4):(i + 5)] > signifs[[gene]]) / 10
#           ) 
#       }
#     }
#     print(ret)
#     ret
#   }
# ) %>% do.call(rbind, .)

gene_signif_freq <- by(wheat_data$snp, wheat_data$snp$chrom,
  function(snp_data) {
    ret <- tibble(
      gene = character(), chrom = integer(), pos_mb = double(), freq = double()
    )
    for (gene in names(signifs)) {
      for (i in 3:(length(snp_data$id) - 2)) {
        num_nearby <- which(
          snp_data$pos_mb >= snp_data$pos_mb[i] - 5 &
          snp_data$pos_mb <= snp_data$pos_mb[i] + 5
        ) %>% length()
        freq <- integer()
        if (num_nearby >= 5) {
          freq <- sum(
            snp_data[[gene]][(i - 2):(i + 2)] > signifs[[gene]]
          ) / 5
        } else {
          freq <- NA
        }
        if (! is.na(freq) && freq == 0) {
          freq <- NA
        }
        ret <- ret %>%
          add_row(
            gene = gene, chrom = snp_data$chrom[1], 
            pos_mb = snp_data$pos_mb[i],
            freq = freq
          ) 
      }
    }
    temp <- vector()
    for (i in 1:nrow(ret)) {
      if (! is.na(ret[i, ]$freq)) {
        temp <- c(temp, i)
      } else if (length(temp) > 0) {
        # highest <- colMeans(ret[temp, "freq"])[[1]]
        highest <- which.max(ret[temp, "freq"][[1]])
        ret[temp[-highest], "freq"] <- NA
        # for (j in temp) {
        #   ret[j, ]$freq <- highest
        # }
        temp <- vector()
      }
    }
    print(ret)
    ret
  }
) %>% do.call(rbind, .)

# load the gene positions
known_genes <- read_csv(
  "Data/Intermediate/Aligned_genes/selected_alignments/known_genes.csv",
  col_names = c("gene", "chrom", "pos", "sleng", "salign", "%id")
) %>%
  select(gene, chrom, pos) %>%
  mutate(pos_mb = pos / 1000000)

# convert known genes chromosome  names to appropriate integer
chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
  t() %>% as.vector()
for (i in seq_along(chroms)) {
  known_genes$chrom[known_genes$chrom == chroms[i]] <- i
}
known_genes$chrom <- as.integer(known_genes$chrom)

# add min phi values of plotting to each gene so they appear at bottom of plots
known_genes <- cbind(known_genes, min_signif = 0)

# add the genes in to the data frame
gene_signif_freq_known <- gene_signif_freq %>%
  rbind.fill(known_genes) %>%
  arrange(chrom, gene, pos) %>%
  as.tibble()

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by Gene or gene type
lables <- c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")
legend_title <- "Comparisons and Genes"
plots <- by(gene_signif_freq_known, gene_signif_freq_known$chrom,
  function(chrom_data) {
    chrom <- chrom_data$chrom[1]
    chrom_data %>%
      ggplot() +
        ylim(0, 1) +
        xlim(0, max_genome_lengths[ifelse(chrom %% 3, chrom %% 3, 3)]) +
        geom_point(
          aes(pos_mb, freq, colour = gene), size = 0.5
        ) +
        geom_text_repel(aes(pos_mb, min_signif, colour = gene,
          label = gene), angle = 90, hjust = 0, vjust = 1, size = 3,
          segment.colour = "black", nudge_y = 0.15, nudge_x = 60,
          show.legend = FALSE
        )
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "Jost D",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Frequency of Nearby Markers in Top 2.5% of Jost D Values for Resistance Gene Groups",
  legend = c(1, 1)
)

plots_matrix[1, 1] <- plots_matrix[1, 1] +
  theme(panel.border = element_rect(fill = NA, colour = "#A3A500", size = 2))
plots_matrix[1, 3] <- plots_matrix[1, 3] +
  theme(panel.border = element_rect(fill = NA, colour = "#00BF7D", size = 2))
plots_matrix[2, 3] <- plots_matrix[2, 3] +
  theme(panel.border = element_rect(fill = NA, colour = "#00B0F6", size = 2))
plots_matrix[5, 3] <- plots_matrix[5, 3] +
  theme(panel.border = element_rect(fill = NA, colour = "#F8766D", size = 2))
plots_matrix[7, 3] <- plots_matrix[7, 3] +
  theme(panel.border = element_rect(fill = NA, colour = "#E76BF3", size = 2))

# plot the matrix
png(str_c("Results/loci/D/known_genes_D_freq_point.png"),
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
print(plots_matrix + theme(legend.position = "bottom"))
dev.off()