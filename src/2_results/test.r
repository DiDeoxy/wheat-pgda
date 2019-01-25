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
extremes <- calc_extremes(wheat_data, groups, prob = 0.95)

# create a table of the regions with a high density of extreme markers
group_extreme_freqs <- calc_group_extreme_freqs(
  wheat_data, extremes
)

comp_gef <- group_extreme_freqs[complete.cases(group_extreme_freqs), ]

all_overlaps <- by(comp_gef, comp_gef$chrom, function (chrom) {
  points <- c(chrom$start, chrom$end) %>% sort()
  # ret <- NULL
  chrom_all_overlaps <- tibble()
  for (i in 1:(length(points) - 1)) {
    # print(paste(points[i], points[i + 1]))
    spannings <- chrom[which(chrom$start <= points[i] & chrom$end >= points[i + 1]), ]
    # print(spannings)
    overlaps_by_group <- list()
    # print(nrow(spannings))
    if (nrow(spannings)) {
      for (j in 1:nrow(spannings)) {
        temp <- tibble(
          chrom = chrom$chrom[1],
          pos_mb = strsplit(spannings[j, "pos_mb"] %>% as.character(), ' ')[[1]] %>% as.numeric(),
          group = spannings$group[j],
          extreme = strsplit(spannings[j, "extreme"] %>% as.character(), ' ')[[1]],
          Ds = strsplit(spannings[j, "Ds"] %>% as.character(), ' ')[[1]],
          ids = strsplit(spannings[j, "ids"] %>% as.character(), ' ')[[1]]
        )
        temp <- temp[which(temp$pos_mb >= points[i] & temp$pos_mb <= points[i + 1]), ]
        names(temp)[3:5] <- str_c(names(temp)[3:5], spannings$group[j], sep = "_")
        overlaps_by_group[[spannings$group[j]]] <- temp
      }
    }
    # print(overlaps_by_group)
    if (length(overlaps_by_group) > 1) {
      chrom_all_overlaps <- bind_rows(chrom_all_overlaps, Reduce((function (x, y) { full_join(x, y) }), overlaps_by_group))
    } else {
      chrom_all_overlaps <- bind_rows(chrom_all_overlaps, overlaps_by_group)
    }
  }
  chrom_all_overlaps
}) %>% Reduce(function (x, y) { bind_rows(x, y) }, .)
all_overlaps <- all_overlaps[,order(colnames(all_overlaps))]

# print our the markers involved in each linked region
genes <- rbind(pheno_genes, resi_genes)
base <- "Results/loci/D/closest_markers"
by(all_overlaps, list(paste(all_overlaps$group_chrs_chrw, all_overlaps$group_chrs_csws, all_overlaps$group_csws_chrw), all_overlaps$chrom), function (chrom_overlap) {
  print(chrom_overlap)
  # file_name <- paste(chrom$chrom[1])
})
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
comp_gef <- group_extreme_freqs[complete.cases(group_extreme_freqs), ]

find_groups <- function(region1, region2) {
  groups <- c(region1[c(2, 3)], region2[c(2, 3)]) %>% unlist() %>% unique()
  if (length(groups) == 1) {
    return (c(groups, NA, NA))
  } else if (length(groups) == 2) {
    return (c(groups, NA))
  } else {
    return (groups)
  }
}

# having problem with double overlaps, insteat of looking for single overlaps first
# look for total overlaps using a which statment follwoing logic of identifying single overlaps
blah <- by(comp_gef, comp_gef$chrom, function (chrom) {
  single_overlaps <- tibble(
    chromo = character(), group1 = character(), group2 = character(),
    start = double(), end = double()
  )
  double_overlaps <- tibble(
    chrom = character(), group1 = character(), group2 = character(),
    group3 = character(), start = double(), end = double()
  )
  if (nrow(chrom) > 1) {
    for (i in 1:(nrow(chrom) - 1)) {
      for (j in (i + 1):nrow(chrom)) {
        if (
          chrom[i, ]$start <= chrom[j, ]$start
          && chrom[i, ]$end >= chrom[j, ]$start
        ) {
          single_overlaps <- single_overlaps %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group,
              group2 = NA, start = chrom[i, ]$start, end = chrom[j, ]$start
            )
          single_overlaps <- single_overlaps %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group,
              group2 = chrom[j, ]$group, start = chrom[j, ]$start,
              end = chrom[i, ]$end
            )
        } else if (
          chrom[i, ]$start >= chrom[j, ]$start
          && chrom[i, ]$end <= chrom[j, ]$end
        ) {
          single_overlaps <- single_overlaps %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[j, ]$group,
              group2 = NA, start = chrom[j, ]$start, end = chrom[i, ]$start
            )
          single_overlaps <- single_overlaps %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group,
              group2 = chrom[j, ]$group, start = chrom[i, ]$start,
              end = chrom[i, ]$end
            )
          single_overlaps <- single_overlaps %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[j, ]$group,
              group2 = NA, start = chrom[i, ]$end, end = chrom[j, ]$end
            )
        } else if (
          chrom[i, ]$start <= chrom[j, ]$end
          && chrom[i, ]$end >= chrom[j, ]$end
        ) {
          single_overlaps <- single_overlaps %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group,
              group2 = chrom[j, ]$group, start = chrom[i, ]$start,
              end = chrom[j, ]$end
            )
          single_overlaps <- single_overlaps %>%
            add_row(
              chromo = chrom$chrom[1], group1 = chrom[i, ]$group, group2 = NA,
              start = chrom[j, ]$end, end = chrom[i, ]$end
            )
        }
        # else {
        #   single_overlaps <- single_overlaps %>%
        #     add_row(
        #       chromo = chrom$chrom[1], group1 = chrom[i, ]$group,
        #       group2 = NA, start = chrom[i, ]$start, end = chrom[i, ]$end
        #     )
        #   break
        # }
      }
    }
    # for (i in 1:(nrow(single_overlaps) - 1)) {
    #   for (j in i:nrow(single_overlaps)) {
    #     groups <- find_groups(single_overlaps[i, ], single_overlaps[j, ])
    #     # print(1)
    #     if (
    #       single_overlaps[i, ]$start <= single_overlaps[j, ]$start
    #       && single_overlaps[i, ]$end >= single_overlaps[j, ]$start
    #     ) {
    #       # print(2)
    #       double_overlaps <- double_overlaps %>%
    #         add_row(
    #           chrom = single_overlaps$chromo[1], group1 = groups[1], group2 = groups[2], 
    #           group3 = groups[3], start = single_overlaps[i, ]$start, end = single_overlaps[j, ]$start
    #         )
    #       double_overlaps <- double_overlaps %>%
    #         add_row(
    #           chrom = single_overlaps$chromo[1], group1 = groups[1], group2 = groups[2], 
    #           group3 = groups[3], start = single_overlaps[j, ]$start,
    #           end = single_overlaps[i, ]$end
    #         )
    #     } else if (
    #       single_overlaps[i, ]$start >= single_overlaps[j, ]$start
    #       && single_overlaps[i, ]$end <= single_overlaps[j, ]$end
    #     ) {
    #       # print(3)
    #       double_overlaps <- double_overlaps %>%
    #         add_row(
    #           chrom = single_overlaps$chromo[1], group1 = groups[1], group2 = groups[2], 
    #           group3 = groups[3], start = single_overlaps[j, ]$start, end = single_overlaps[i, ]$start
    #         )
    #       double_overlaps <- double_overlaps %>%
    #         add_row(
    #           chrom = single_overlaps$chromo[1], group1 = groups[1], group2 = groups[2], 
    #           group3 = groups[3], start = single_overlaps[i, ]$start, end = single_overlaps[i, ]$end
    #         )
    #       double_overlaps <- double_overlaps %>%
    #         add_row(
    #           chrom = single_overlaps$chromo[1], group1 = groups[1], group2 = groups[2], 
    #           group3 = groups[3], start = single_overlaps[i, ]$end, end = single_overlaps[j, ]$end
    #         )
    #     } else if (
    #       single_overlaps[i, ]$start <= single_overlaps[j, ]$end
    #       && single_overlaps[i, ]$end >= single_overlaps[j, ]$end
    #     ) {
    #       # print(4)
    #       double_overlaps <- double_overlaps %>%
    #         add_row(
    #           chrom = single_overlaps$chromo[1], group1 = groups[1], group2 = groups[2], 
    #           group3 = groups[3], start = single_overlaps[i, ]$start, end = single_overlaps[j, ]$end
    #         )
    #       double_overlaps <- double_overlaps %>%
    #         add_row(
    #           chrom = single_overlaps$chromo[1], group1 = groups[1], group2 = groups[2], 
    #           group3 = groups[3], start = single_overlaps[j, ]$end, end = single_overlaps[i, ]$end
    #         )
    #     }
    #     # else {
    #     #   # print(5)
    #     #   groups <- find_groups(single_overlaps[i, ], single_overlaps[i, ])
    #     #   double_overlaps <- double_overlaps %>%
    #     #     add_row(
    #     #       chrom = single_overlaps$chromo[1], group1 = groups[1], group2 = groups[2], 
    #     #       group3 = groups[3], start = single_overlaps[i, ]$start, end = single_overlaps[i, ]$end
    #     #     )
    #     # }
    #   }
    # }
  }
  single_overlaps
  # double_overlaps
})
blah

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
# install.packages("data.table")

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

comp_gef <- group_extreme_freqs[complete.cases(group_extreme_freqs), ]
blah <- by(comp_gef, comp_gef$chrom, function (chrom) {
  chrom <- chrom %>% arrange(start)
  for (i in 1:nrow(chrom)) {
    ovlp <- which(chrom$start <= chrom$end[i])
    if (any(ovlp > i)) {
      for (j in 1:nrow(chrom[ovlp]))
      }
    }
  }
})