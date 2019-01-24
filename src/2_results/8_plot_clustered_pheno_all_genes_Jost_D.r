library(plyr)
library(tidyverse)
library(GGally)
library(ggrepel)
library(extrafont)

# load custom functions
source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/colour_sets.R")
source("src/R_functions/funcs_locus_by_locus.R")

# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

# find the max position of any marker on each genome for xlims
max_genome_lengths <- calc_max_genome_lengths(wheat_data)

# names of groups to be plotted
groups <- c("chrs_csws", "chrs_chrw", "csws_chrw")

# find the jost's D values of each marker in each Gene and add to data set
wheat_data <- add_group_stat(wheat_data, groups)

# find the extreme threshold for each Gene
extremes <- calc_extremes(wheat_data, groups, prob = 0.955)

# create a table of the regions with a high density of extreme markers
group_extreme_freqs <- calc_group_extreme_freqs(
  wheat_data, extremes
)

# load the gene positions
pheno_genes <- load_groups("pheno_genes.csv", base = 0.25) %>%
  mutate(group = "pheno_gene") %>%
  dplyr::rename(mean_pos_mb = "pos_mb")
resi_genes <- load_groups("resi_genes.csv", base = 0.25) %>%
  mutate(group = "resi_gene") %>%
  dplyr::rename(mean_pos_mb = pos_mb)

# add the genes positons to the regions table
group_extreme_freqs_genes <- group_extreme_freqs %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  arrange(chrom, group, pos_mb) %>%
  as.tibble()

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
lables <- c(
  "CHRS vs CHRW", "CHRS vs CSWS", "CSWS vs CHRW", "Phenotype Genes",
  "Resistance Genes"
)
legend_title <- "Comparisons and Genes"
group_extreme_freqs_genes$mean_pos_mb[is.na(which(group_extreme_freqs_genes$base))]
plots <- by(
  group_extreme_freqs_genes, group_extreme_freqs_genes$chrom,
  function(chrom_data) {
    chrom <- chrom_data$chrom[1]
    chrom_data %>%
      ggplot() +
      ylim(0.25, 1) +
      xlim(
        0, 
        max_genome_lengths[
          ifelse(grepl("A", chrom), 1, ifelse(grepl("B", chrom), 2, 3))
        ]
      ) +
      # geom_segment(
      #   aes(
      #     x = start - 1, y = mean_D, xend = end + 1, yend = mean_D,
      #     size = num_linked
      #   )
      # ) +
      # scale_size_continuous(
      #   "Number of Markers", trans = "sqrt",
      #   limits = c(1, max(group_extreme_freqs$num_linked, na.rm = T)),
      #   breaks = c(25, 50, 100, 200)
      # ) +
      geom_segment(
        aes(
          x = start - 1, y = mean_D, xend = end + 1, yend = mean_D,
          colour = group
        ), size = 1
      ) +
      geom_point(
        aes(start, mean_D, colour = group), shape = 20, size = 2
      ) +
      geom_point(
        aes(end, mean_D, colour = group), shape = 20, size = 2
      ) +
      geom_point(
        aes(mean_pos_mb, base, colour = group, shape = group), size = 1
      ) +
      geom_text_repel(
        aes(mean_pos_mb, base, colour = group, label = id), angle = 90, hjust = 0,
        vjust = 1, size = 3, segment.colour = "black",
        nudge_y = 0.07,
        nudge_x = ifelse(chrom == "1D", 80,
          ifelse(chrom %in% c("2D", "4A"), -60, 40)
        ),
        show.legend = FALSE
      ) +
      scale_colour_manual(
        legend_title, labels = lables, values = colours_groups_genes,
        limits = levels(as.factor(group_extreme_freqs_genes$group))
      ) +
      scale_shape_manual(
        legend_title, labels = lables, values = c(15, 16, 18, 17, 8),
        limits = levels(as.factor(group_extreme_freqs_genes$group))
      ) +
      labs(colour = "Comparison")
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = str_c(
    "Average Jost's D Values in Regions with Nearby Exceptional Markers"
  ),
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Number of Markers in a Region and their Average Jost's D Values",
  legend = c(1, 1)
)

# plot the matrix
png(str_c("Results/loci/D/comps_D.png"),
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()

# # print our the markers involved in each linked region
comp_gef <- group_extreme_freqs[complete.cases(group_extreme_freqs), ]
genes <- rbind(pheno_genes, resi_genes)
base <- "Results/loci/D/closest_markers"
for (row in 1:nrow(comp_gef)) {
  file_name1 <- paste(
    cbind(
      comp_gef[row, c("chrom", "mean_pos_mb", "group", "num_linked")] %>%
      round_df(0),
      comp_gef[row, c("freq_extreme", "mean_D")] %>% round_df(2)
    ), collapse = '_'
  )
  file_name2 <- paste(
    cbind(
      comp_gef[row, c("group", "chrom", "mean_pos_mb", "num_linked")] %>%
      round_df(0),
      comp_gef[row, c("freq_extreme", "mean_D")] %>% round_df(2)
    ), collapse = '_'
  )
  linked <- tibble(
    extreme = strsplit(comp_gef[row, "extreme"] %>% as.character(), ' ')[[1]],
    pos_mb = strsplit(comp_gef[row, "pos_mb"] %>% as.character(), ' ')[[1]] %>% as.numeric(),
    Ds = strsplit(comp_gef[row, "Ds"] %>% as.character(), ' ')[[1]],
    ids = strsplit(comp_gef[row, "ids"] %>% as.character(), ' ')[[1]]
  )
  for (row2 in 1:nrow(genes)) {
    if (genes[row2, ]$chrom == comp_gef[row, ]$chrom) {
      linked <- linked %>%
        add_row(
          extreme = "NA", pos_mb = genes[row2, ]$mean_pos_mb, Ds = "NA",
          ids = genes[row2, ]$id
        )
    }
  }
  linked <- linked %>% arrange(pos_mb)
  ifelse(! dir.exists(file.path(base)), dir.create(file.path(base)), FALSE)
  write_csv(linked, file.path(base, str_c(file_name1, ".csv")))
  write_csv(linked, file.path(base, str_c(file_name2, ".csv")))
}

################################################################################
# source("http://bioconductor.org/biocLite.R")
# biocLite("IRanges")
# library(IRanges)

# comp_gef <- group_extreme_freqs[complete.cases(group_extreme_freqs), ]
# ?IRanges
# comp_gef
# by(comp_gef, comp_gef$chrom, function (chrom) {
#   ir <- IRanges(
#     start = chrom$start, end = chrom$end,
#     # names = comp_gef$group
#   )
# })
# ir <- IRanges(
#   start = comp_gef$start, end = comp_gef$end,
#   # names = str_c(comp_gef$chrom, comp_gef$group, collapse = " ")
# )
install.packages("data.table")

# num_chrs_chrw_chrs_csws <- 0
# num_chrs_chrw_csws_chrw <- 0
# num_chrs_csws_csws_chrw <- 0
# blah <- by(comp_gef, comp_gef$chrom, function (chrom) {
#   if (nrow(chrom) > 1) {
#     pairs <- list()
#     for (i in 1:(nrow(chrom) - 1)) {
#       for (j in (i + 1):nrow(chrom)) {
#         if (1 %in%
#           findInterval(
#             c(chrom[i, ]$start, chrom[i, ]$end),
#             c(chrom[j, ]$start, chrom[j, ]$end)
#           )
#           ||
#           (
#             c(chrom[i, ]$start, chrom[i, ]$end) ==
#             c(chrom[j, ]$start, chrom[j, ]$end)
#           )
#         ) {
#           if (
#             (chrom[i, ]$group == "chrs_chrw" && chrom[j, ]$group == "chrs_csws")
#             ||
#             (chrom[i, ]$group == "chrs_csws" && chrom[j, ]$group == "chrs_chrw")
#           ) {
#             num_chrs_chrw_chrs_csws <<- num_chrs_chrw_chrs_csws + 1
#             pairs[["chrs_chrw_chrs_csws"]][length(pairs[["chrs_chrw_chrs_csws"]]) + 1] <- paste(chrom$chrom[1], mean(chrom[i, ]$start, chrom[i, ]$end), collapse = " ")
#           } else if (
#             (chrom[i, ]$group == "chrs_chrw" && chrom[j, ]$group == "csws_chrw")
#             ||
#             (chrom[i, ]$group == "csws_chrw" && chrom[j, ]$group == "chrs_chrw")
#           ) {
#             num_chrs_chrw_csws_chrw <<- num_chrs_chrw_csws_chrw + 1
#             pairs[["chrs_chrw_csws_chrw"]][length(pairs[["chrs_chrw_csws_chrw"]]) + 1] <- paste(chrom$chrom[1], mean(chrom[i, ]$start, chrom[i, ]$end), collapse = " ")
#           } else {
#             num_chrs_csws_csws_chrw <<- num_chrs_csws_csws_chrw + 1
#             pairs[["chrs_csws_csws_chrw"]][length(pairs[["chrs_csws_csws_chrw"]]) + 1] <- paste(chrom$chrom[1], mean(chrom[i, ]$start, chrom[i, ]$end), collapse = " ")
#           }
#         }
#       }
#     }
#   }
#   pairs
# })
# num_chrs_chrw_chrs_csws
# num_chrs_chrw_csws_chrw
# num_chrs_csws_csws_chrw
# blah

# test <- function(regions, row) {
#   (
#     sum(regions[row - 1, 3:5] %in% regions[row, 3:5]) == 3
#     && regions[row - 1, ]$interval == regions[row, ]$interval - 1
#   )
# }

# comp_gef <- group_extreme_freqs[complete.cases(group_extreme_freqs), ]
# comp_gef %>% print(n = Inf)
# blah <- by(comp_gef, comp_gef$chrom, function (chrom) {
#   regions <- tibble(
#     chrom = character(), interval = integer(), group1 = character(),
#     group2 = character(), group3 = character()
#   )
#   if (nrow(chrom) > 1) {
#     for (i in 1:max(chrom$end %>% as.integer() + 1)) {
#       groups <- NULL
#       for (row in 1:nrow(chrom)) {
#         if (chrom_start >= i - 1
#           chrom[row, ]$start <= i  && i <= chrom[row, ]$end) {
#           groups <- c(groups, chrom[row, ]$group)
#         }
#       }
#       if (length(groups)) {
#         if (length(groups) == 1) {
#           regions <- regions %>%
#             add_row(
#               chrom = chrom$chrom[1], interval = i, group1 = groups[1]
#             )
#         } else if (length(groups) == 2) {
#           regions <- regions %>%
#             add_row(
#               chrom = chrom$chrom[1], interval = i, group1 = groups[1],
#               group2 = groups[2]
#             )
#         } else {
#           regions <- regions %>%
#             add_row(
#               chrom = chrom$chrom[1], interval = i, group1 = groups[1],
#               group2 = groups[2],group3 = groups[3]
#             )
#         }
#       }
#     }
#   }
#   over_ints <- tibble(
#     chrom = character(), start = integer(), end = integer(),
#     group1 = character(), group2 = character(), group3 = character()
#   )
#   if (nrow(regions)) {
#     start <- regions[1, ]$interval
#     end <- regions[1, ]$interval
#     for (row in 2:nrow(regions)) {
#       if (test(regions, row)) {
#         end <- regions[row, ]$interval
#       }
#       if (row == nrow(regions) || ! test(regions, row)) {
#         over_ints <- over_ints %>%
#           add_row(
#             chrom = chrom$chrom[1], start = start, end = end,
#             group1 = regions[row - 1, ]$group1,
#             group2 = regions[row - 1, ]$group2,
#             group3 = regions[row - 1, ]$group3
#           )
#         start <- regions[row, ]$interval
#         end <- regions[row, ]$interval
#       }
#     }
#   }
#   over_ints
# })
# for (ele in 1:length(blah)) {
#   blah[[ele]] %>% print(n = Inf)
# }
# blah[[length(blah)]] %>% print(n = Inf)
# num_chrs_chrw_chrs_csws
# num_chrs_chrw_csws_chrw
# num_chrs_csws_csws_chrw

find_groups <- function(region1, region2) {
  groups <- c(region1[c(2, 3)], region2[c(2, 3)]) %>% unique()
  if (length(groups) == 1) {
    return (c(groups, NA, NA))
  } else if (length(groups) == 2) {
    return (c(groups, NA))
  } else {
    return (groups)
  }
}

blah <- by(comp_gef, comp_gef$chrom, function (chrom) {
  overlaps_two <- tibble(
    chromo = character(), group1 = character(), group2 = character(),
    start = double(), end = double()
  )
  if (nrow(chrom) > 1) {
    for (i in 1:(nrow(chrom) - 1)) {
      for (j in (i + 1):nrow(chrom)) {
        if (
          chrom[i, ]$start <= chrom[j, ]$start
          && chrom[i, ]$end >= chrom[j, ]$start
        ) {
          overlaps_two <- overlaps_two %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group,
              group2 = NA, start = chrom[i, ]$start, end = chrom[j, ]$start
            )
          overlaps_two <- overlaps_two %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group,
              group2 = chrom[j, ]$group, start = chrom[j, ]$start,
              end = chrom[i, ]$end
            )
        } else if (
          chrom[i, ]$start >= chrom[j, ]$start
          && chrom[i, ]$end <= chrom[j, ]$end
        ) {
          overlaps_two <- overlaps_two %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[j, ]$group,
              group2 = NA, start = chrom[j, ]$start, end = chrom[i, ]$start
            )
          overlaps_two <- overlaps_two %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group,
              group2 = chrom[j, ]$group, start = chrom[i, ]$start,
              end = chrom[i, ]$end
            )
          overlaps_two <- overlaps_two %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[j, ]$group,
              group2 = NA, start = chrom[i, ]$end, end = chrom[j, ]$end
            )
        } else if (
          chrom[i, ]$start <= chrom[j, ]$end
          && chrom[i, ]$end >= chrom[j, ]$end
        ) {
          overlaps_two <- overlaps_two %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group,
              group2 = chrom[j, ]$group, start = chrom[i, ]$start,
              end = chrom[j, ]$end
            )
          overlaps_two <- overlaps_two %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group, group2 = NA,
              start = chrom[j, ]$end, end = chrom[i, ]$end
            )
        }
      }
    }
    overlaps_three <- tibble(
      chrom = character(), group1 = character(), group2 = character(),
      group3 = character(), start = double(), end = double()
    )
    for (i in 1:(nrow(overlaps_two) - 1)) {
      for (j in i:nrow(overlaps_two)) {
        groups <- find_groups(overlaps_two[i, ], overlaps_two[j, ])
        if (
          single_overlaps[i, ]$start <= overlaps_two[j, ]$start
          && overlaps_two[i, ]$end >= overlaps_two[j, ]$start
        ) {
          overlaps_three <- overlaps_three %>%
            add_row(
              chrom = overlaps_two$chromo[1], group1 = groups[1], group2 = groups[2], 
              group3 = groups[3], start = overlaps_two[i, ]$start, end = overlaps_two[j, ]$start
            )
          overlaps_three <- overlaps_three %>%
            add_row(
              chrom = overlaps_two$chromo[1], group1 = groups[1], group2 = groups[2], 
              group3 = groups[3], start = overlaps_two[j, ]$start,
              end = overlaps_two[i, ]$end
            )
        } else if (
          overlaps_two[i, ]$start >= overlaps_two[j, ]$start
          && overlaps_two[i, ]$end <= overlaps_two[j, ]$end
        ) {
          overlaps_three <- overlaps_three %>%
            add_row(
              chrom = overlaps_two$chromo[1], group1 = groups[1], group2 = groups[2], 
              group3 = groups[3], start = overlaps_two[j, ]$start, end = overlaps_two[i, ]$start
            )
          overlaps_three <- overlaps_three %>%
            add_row(
              chrom = overlaps_two$chromo[1], group1 = groups[1], group2 = groups[2], 
              group3 = groups[3], start = overlaps_two[i, ]$start, end = overlaps_two[i, ]$end
            )
          overlaps_three <- overlaps_three %>%
            add_row(
              chrom = overlaps_two$chromo[1], group1 = groups[1], group2 = groups[2], 
              group3 = groups[3], start = overlaps_two[i, ]$end, end = overlaps_two[j, ]$end
            )
        } else if (
          overlaps_two[i, ]$start <= overlaps_two[j, ]$end
          && overlaps_two[i, ]$end >= overlaps_two[j, ]$end
        ) {
          overlaps_three <- overlaps_three %>%
            add_row(
              chrom = overlaps_two$chromo[1], group1 = groups[1], group2 = groups[2], 
              group3 = groups[3], start = overlaps_two[i, ]$start, end = overlaps_two[j, ]$end
            )
          overlaps_three <- overlaps_three %>%
            add_row(
              chrom = overlaps_two$chromo[1], group1 = groups[1], group2 = groups[2], 
              group3 = groups[3], start = overlaps_two[j, ]$end, end = overlaps_two[i, ]$end
            )
        }
      }
    }
  }
  overlaps_three
})
blah

sum(comp_gef$group == "chrs_chrw")
comp_gef[which(comp_gef$group == "chrs_chrw"), ] %>% print(n = Inf)



library(adegenet)
comp_genind <- read_rds("Data/Intermediate/Adegenet/chrs_chrw_genind.rds")
markers <- c(least = "BS00071087_51", most = "wsnp_Ex_c41347_48190370")
for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind@pop, allele = comp_genind@tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}
comp_genind <- read_rds("Data/Intermediate/Adegenet/chrs_csws_genind.rds")
markers <- c("Kukri_c10033_724")
for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind@pop, allele = comp_genind@tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}
# sum(comp_gef$group == "chrs_chrw" & comp_gef$mean_D > 0.5)
# sum(comp_gef$group == "chrs_chrw" & comp_gef$mean_D > 0.75)
comp_gef[which(comp_gef$group == "chrs_chrw" & comp_gef$mean_D > 0.75), ] %>% print(n = Inf)
comp_gef[which(comp_gef$group == "chrs_chrw" & comp_gef$num_linked > 20), ] %>% print(n = Inf)
# comp_gef[comp_gef$group == "chrs_chrw", ] %>% print(n = Inf)
sum(comp_gef$group == "chrs_csws")
sum(comp_gef$group == "csws_chrw")


gef_chrs_chrw <- comp_gef[which(comp_gef$group == "chrs_chrw"), ]
total_extreme <- 0
total_D <- 0
for (i in 1:nrow(gef_chrs_chrw)) {
  total_extreme <- total_extreme + gef_chrs_chrw[i, ]$num_linked * gef_chrs_chrw[i, ]$freq_extreme
  total_D <- total_D + gef_chrs_chrw[i, ]$num_linked * gef_chrs_chrw[i, ]$mean_D
}
total_extreme/sum(gef_chrs_chrw$num_linked)
total_D/sum(gef_chrs_chrw$num_linked)