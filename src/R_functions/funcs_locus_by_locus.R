library(tidyverse)

calc_eh <- function (genotypes) {
  apply(genotypes, 1, function (snp) {
    2 * ((sum(snp == 0) / sum(snp == 0 | 2)) *
      (sum(snp == 2) / sum(snp == 0 | 2)))
  })
}

calc_max_genome_lengths <- function(wheat_data) {
  by(wheat_data$snp$pos_mb, wheat_data$snp$chrom, max) %>% 
    (
      function (max_chrom_lengths) {
        c(
          A = max(max_chrom_lengths[seq(1, 19, 3)]),
          B = max(max_chrom_lengths[seq(2, 20, 3)]),
          D = max(max_chrom_lengths[seq(3, 21, 3)])
        )
      }
    )
}

calc_max_homeolog_lengths <- function(wheat_data) {
  by(wheat_data$snp$pos, wheat_data$snp$chrom, max) %>%
    (
      function (max_chrom_lengths) {
        c(
          one = max(max_chrom_lengths[c(1, 2, 3)]),
          two = max(max_chrom_lengths[c(4, 5, 6)]),
          three = max(max_chrom_lengths[c(7, 8, 9)]),
          four = max(max_chrom_lengths[c(10, 11, 12)]),
          five = max(max_chrom_lengths[c(13, 14, 15)]),
          six = max(max_chrom_lengths[c(16, 17, 18)]),
          seven = max(max_chrom_lengths[c(19, 20, 21)])
        )
      }
    )
}

add_group_stat <- function(wheat_data, groups) {
  for (group in groups) {
    wheat_data$snp <- wheat_data$snp %>% 
      add_column(
        !!group := read_rds(
          str_c("Data/Intermediate/mmod/", group,"_Jost_D.rds")
        )[[1]]
      )
  }
  wheat_data
}
group <- "Lr34"

calc_extremes <- function(wheat_data, groups, prob = 0.975) {
  temp <- vector()
  for (group in groups) {
    temp[group] <- quantile(wheat_data$snp[[group]], prob = prob, na.rm = TRUE)
  }
  temp
}

load_groups <- function(csv, base = 0) {
  read_csv(
    str_c("Data/Intermediate/Aligned_genes/selected_alignments/", csv),
    col_names = c("id", "chrom", "pos", "sleng", "salign", "%id")
  ) %>%
    mutate(pos_mb = pos / 1e6) %>%
    select(id, chrom, pos_mb) %>%
    cbind(base = base)
}

snp_densities <- function(chrom) {
  # the distances between markers on the chromosome
  gaps <- diff(chrom$pos_mb)
  # a vector for holding the average density of snps near snp i
  densities <- vector()
  for (i in 1:nrow(chrom)) {
    # a vector for holding the inices of the 10 snps upstream and 10 snps
    # downstream from snp i, as well as snp i
    nearest <- vector()
    if (i < 5 && length(gaps) >= 8) {
      # position 1 in gaps is the gap between snps 1 and 2
      # position 20 in gaps is the gap between snps 20 and 21
      nearest <- 1:(i + 4)
    } else if (i > length(gaps) - 5 && length(gaps) >= 8) {
      # position length(gaps) - 20 in gaps is the gap between snps
      # length(chrom) - 21  and length(chrom) - 20
      # position length(gaps) in gaps is the gap between snp length(chrom) - 1
      # and length(chrom)
      nearest <- (i - 4):length(gaps)
    } else if (length(gaps) < 8) {
      nearest <- 1:length(gaps)
    } else {
      # the positon i - 9 in gaps is the gap between snps i - 10 and i - 9
      # the psotion i + 9 in gaps is the gap between snps i + 9 and i + 10
      nearest <- (i - 4):(i + 4)
    }
    # calc the density of the regions aorund marker i and add it to the
    # densities vector
    densities <- c(densities, mean(gaps[nearest]))
  }
  densities
}

find_windows <- function (snp_data, extremes) {
  # initialize an empty tibble with the needed columns
  windows <- tibble(
    chrom = character(), group = character(), pos_mb = double(),
    freq_nearby_extreme = double(), D = double(), id = character()
  )
  # the max distance in Mb from the considered marker that the upstream and
  # downstream markers of num can be
  dist <- 15
  # for each snp in the comparison
  for (i in 1:nrow(snp_data)) {
    # for holding those markers considered nearby (within num & dist)
    nearby <- i
    # find upstream nearby markers
    if (i > 1) {
      for (j in (i - 1):1) {
        if (
          j >= 1
          && snp_data$pos_mb[i] - snp_data$pos_mb[j] <= dist
          # && i - j <= num
        ) {
          nearby <- c(j, nearby)
        } else {
          break
        }
      }
    }
    # find downstream nearby markers
    if (i < nrow(snp_data)) {
      for (j in (i + 1):nrow(snp_data)) {
        if (
          j <= nrow(snp_data)
          && snp_data$pos_mb[j] - snp_data$pos_mb[i] <= dist
          # && j - i <= num
        ) {
          nearby <- c(nearby, j)
        } else {
          break
        }
      }
    }
    # for each comparison
    for (group in names(extremes)) {
      # store the useful data for each marker in each group including its
      # frequency of nearby extreme markers, its Jost's D values, and its id
      windows <- windows %>%
        add_row(
          chrom = snp_data$chrom[i], group = group, pos_mb = snp_data$pos_mb[i],
          freq_nearby_extreme = sum(
            snp_data[[group]][nearby] > extremes[group]
          ) / length(nearby),
          D = snp_data[[group]][i], id = snp_data$id[i]
        )
    }
  }
  windows
}

find_regions <- function (windows, extremes, freq) {
  linked <- vector()
  regions <- tibble(
    chrom = character(), start = double(), end = double(), group = character(),
    num_linked = integer(), num_extreme = double(), mean_D = double(),
    extreme = character(), pos_mb = character(), Ds = character(),
    ids = character()
  )
  by(windows, list(windows$group, windows$chrom), function (group_chrom) {
    # print(group_chrom)
    # makes sure each group is represented on each chromosome
    group_chrom_regions <- regions %>% add_row(chrom = group_chrom$chrom[1], group = group_chrom$group[1])
    for (i in 1:nrow(group_chrom)) {
      # if a marker has nearby extreme markers it is added to the linked vector
      if (group_chrom$freq_nearby_extreme[i] >= freq) {
        linked <- c(linked, i)
      }
      # needs to be a separate else because the last row can be in linked
      # print(c(group_chrom$freq_nearby_extreme[i], ))
      if (
        (
          ! group_chrom$freq_nearby_extreme[i]
          || i == nrow(group_chrom)
        )
        && length(linked)
        && any(group_chrom$D[linked] > extremes[group_chrom$group[i]])
      ) {
        linked_extreme <- which(group_chrom$D[linked] > extremes[group_chrom$group[i]])
        linked_pruned <- linked[
          linked_extreme[1]:linked_extreme[length(linked_extreme)]
        ]
        linked_pruned <- linked[
          linked_extreme[1]:linked_extreme[length(linked_extreme)]
        ]
        region_extreme <- group_chrom$D[linked_pruned] > extremes[group_chrom$group[i]]
        group_chrom_regions <- group_chrom_regions %>% add_row(
          chrom = group_chrom$chrom[i], group = group_chrom$group[i],
          start = group_chrom$pos_mb[linked_pruned[1]],
          end = 
            group_chrom$pos_mb[linked_pruned[length(linked_pruned)]],
          num_linked = length(linked_pruned),
          num_extreme = sum(region_extreme ),
          mean_D = mean(group_chrom$D[linked_pruned]),
          extreme = paste(region_extreme, collapse = ' '),
          pos_mb = paste(
            group_chrom$pos_mb[linked_pruned], collapse = ' '
          ),
          Ds = paste(group_chrom$D[linked_pruned], collapse = ' '),
          ids = paste(group_chrom$id[linked_pruned], collapse = ' ')
        )
        linked <- vector()
      }
    }
    group_chrom_regions %>% print(n = Inf)
    group_chrom_regions
  }) %>% do.call(rbind, .)
}

calc_group_extreme_freqs <- function(wheat_data, extremes, freq) {
  by(wheat_data$snp, wheat_data$snp$chrom,
    function(snp_data) {
      find_regions(find_windows(snp_data, extremes), extremes, freq)
    }
  ) %>% do.call(rbind, .)  
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits) * 10^digits

  (df)
}


all_ovlps_markers <- function (comp_gef) {
  ovlp_grps <- list()
  all_overlaps <- by(comp_gef, comp_gef$chrom, function (chrom) {
    # points <- c(chrom$start, chrom$end) %>% sort()
    points <- c(chrom$start, chrom$end) %>% unique() %>% sort()
    chrom_all_overlaps <- tibble()
    for (i in 1:(length(points) - 1)) {
      spannings <- chrom[which(chrom$start <= points[i] & chrom$end >= points[i + 1]), ]
      if (! nrow(spannings)) {
        spannings <- chrom[which(chrom$start == points[i] & chrom$end == points[i]), ]
      }
      ovlp_grps[[length(ovlp_grps) + 1]] <<- spannings$group
      overlaps_by_group <- list()
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
      if (length(overlaps_by_group) > 1) {
        chrom_all_overlaps <- bind_rows(chrom_all_overlaps, Reduce((function (x, y) { full_join(x, y) }), overlaps_by_group))
      } else {
        chrom_all_overlaps <- bind_rows(chrom_all_overlaps, overlaps_by_group)
      }
    }
    chrom_all_overlaps
  }) %>% Reduce(function (x, y) { bind_rows(x, y) }, .)
  all_overlaps[, order(colnames(all_overlaps))][, c(1, 12, 5:7, 8:10, 2:4, 11)]
}

print_ovlps_by_chrom <- function(all_overlaps) {
  # print our the markers involved in each linked region
  genes <- rbind(pheno_genes, resi_genes) %>% as.tibble()
  base <- "Results/loci/D/closest_markers"
  ifelse(! dir.exists(file.path(base)), dir.create(file.path(base)), FALSE)
  blah <- by(all_overlaps, all_overlaps$chrom, function (chrom_overlaps) {
    for (row in 1:nrow(genes)) {
      if (genes[row, ]$chrom == chrom_overlaps$chrom[1]) {
        chrom_overlaps <- chrom_overlaps %>%
          add_row(
            chrom = genes[row, ]$chrom, pos_mb = genes[row, ]$mean_pos_mb, ids = genes[row, ]$id
          )
      }
    }
    chrom_overlaps <- chrom_overlaps %>% dplyr::arrange(pos_mb)
    file_name <- paste0(chrom_overlaps$chrom[1], "_overlaps_genes_details.csv")
    file_conn <- file(file.path(base, file_name))
    writeLines(
      c(
        paste(names(chrom_overlaps), collapse = ","),
        paste(chrom_overlaps[1, ], collapse = ",")
      ), file_conn
    )
    close(file_conn)
    file_conn <- file(file.path(base, file_name), open = "at")
    for (row in 2:nrow(chrom_overlaps)) {
      test1 <- chrom_overlaps[row - 1, 6:8]
      test1[is.na(test1)] <- "None"
      test2 <- chrom_overlaps[row , 6:8]
      test2[is.na(test2)] <- "None"
      if (
        ! all(test1 == test2)
      ) {
        writeLines("", file_conn)
        writeLines(paste(chrom_overlaps[row, ], collapse = ","), file_conn)
      } else if (
        chrom_overlaps[row, ]$pos_mb - chrom_overlaps[row - 1, ]$pos_mb > 30
      ) {
        writeLines(rep("", 3), file_conn)
        writeLines(paste(chrom_overlaps[row, ], collapse = ","), file_conn)
      } else {
        writeLines(paste(chrom_overlaps[row, ], collapse = ","), file_conn)
      }
    }
    close(file_conn)
  })
}

comp_ovlps <- function (comp_gef, comp) {
  by(comp_gef, comp_gef$chrom, function (chrom) {
    all_ovlps <- list()
    group_index <- which(chrom$group == comp)
    if (nrow(chrom) > 1 && length(group_index) && length(group_index) < nrow(chrom)) {
      group <- chrom[group_index, ]
      not_group <- chrom[-group_index, ]
      for (i in 1:nrow(group)) {
        ovlps <- group$group[i]
        for (j in 1:nrow(not_group)) {
          if (
            (
              group$start[i] <= not_group$start[j] 
              && (not_group$end[j] <= group$end[i] || not_group$start[j] <= group$end[i])
            )
            ||
            (
              not_group$start[j] <= group$start[i] 
              && (group$end[i] <= not_group$end[j] || group$start[i] <= not_group$end[j])
            )
          ) {
            ovlps <- c(ovlps, not_group$group[j])
          }
        }
        all_ovlps[[i]] <- sort(ovlps) %>% paste(collapse = " ")
      }
    } else {
      for (i in 1:nrow(chrom)) {
        if (chrom$group[i] == comp) {
          all_ovlps[[i]] <- chrom$group[i]
        }
      }
    }
    all_ovlps
  })
}