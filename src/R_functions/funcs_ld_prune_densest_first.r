library(SNPRelate)
library(igraph)
library(tidyverse)

load_data_internal <- function (wheat, maf) {
  print(2)
  # extract the snp metadata for the selected snp
  snp_id <- as.character(
    read.gdsn(index.gdsn(wheat, "snp.id")))
  snp_chrom <- as.integer(
    read.gdsn(index.gdsn(wheat, "snp.chromosome")))
  snp_pos <- as.integer(
    read.gdsn(index.gdsn(wheat, "snp.position"))) / 1000000
  # find the subset of snps with maf above the given threshold
  snp_maf <- snpgdsSelectSNP(wheat, autosome.only = F, maf = maf)
  snp_maf_indices <- match(unlist(snp_maf), snp_id)
  # create a tibble of the snp metadata
  tibble(
    id = snp_id[snp_maf_indices], chrom = snp_chrom[snp_maf_indices],
    pos = snp_pos[snp_maf_indices])
}

calc_genome_mean_dist <- function (snp_data) {
  print(3)
  mean(unlist(by(snp_data, snp_data$chrom, function (chrom) {
    diff(chrom$pos)
  })))
}

calc_coverage <- function (regions) {
  coverage <- 0
  start <- 0
  end <- 0
  for (region in 1:nrows(regions)) {
    if (! start) {
      start <- region$start
      end <- region$end
      next
    }
    if (region$start < end) {
      end <- region$start
    } else {
      coverage <- coverage + end - start
      start <- 0
      end <- 0
    }
  }
  coverage
}

find_clique_regions <- function (centers, pos) {
  centers <- sort(centers)
  range_between_centers <- mean(diff(centers)) / 2
  regions <- tibble()
  for (center in centers) {
    if (center - range_between_centers > 1 &&
      center + range_between_centers < pos[length(pos)]) {
        regions <- rbind(regions,
          c(center - range_between_centers, center + range_between_centers)
        )
      } else if (center - range_between_centers > 1) {
        regions <- rbind(regions,
          c(center - range_between_centers, pos[length(pos)]))
      } else {
        regions <- rbind(regions, c(1, center + range_between_centers))
      }
  }
  colnames(regions) <- c(start, end)
  regions
}

snp_densities <- function(chrom) {
  print(4)
  # the distances between markers on the chromosome
  gaps <- diff(chrom$pos)
  # a vector for holding the average density of snps near snp i
  densities <- vector()
  for (i in 1:length(gaps)) {
    # a vector for holding the inices of the 10 snps upstream and 10 snps
    # downstream from snp i, as well as snp i
    nearest <- vector()
    if (i < 10 && length(gaps) >= 20) {
      # print("a")
      # position 1 in gaps is the gap between snps 1 and 2
      # position 20 in gaps is the gap between snps 20 and 21
      nearest <- 1:20
    } else if (i > length(gaps) - 10 && length(gaps) >= 20) {
      # print("b")
      # position length(gaps) - 20 in gaps is the gap between snps
      # length(chrom) - 21  and length(chrom) - 20
      # position length(gaps) in gaps is the gap between snp length(chrom) - 1
      # and length(chrom)
      nearest <- (length(gaps) - 20):length(gaps)
      # print(c("i: ", i, " nearest: ", nearest))
    } else if (length(gaps) < 20) {
      nearest <- 1:length(gaps)
    } else {
      # print("c")
      # the positon i - 9 in gaps is the gap between snps i - 10 and i - 9
      # the psotion i + 9 in gaps is the gap between snps i + 9 and i + 10
      nearest <- (i - 9):(i + 9)
    }
    # calc the density of the regions aorund marker i and add it to the
    # densities vector
    densities <- c(densities, mean(gaps[nearest]))
  }

  densities
}

largest_clique_finder <- function(chrom, sub_mat, ld_floor, ld_ceiling) {
  print(5)
  # a function used by mapply to calcualte the adjacency of each pair of snps
  # in sub_mat based on proximity and ld
  adj_vec_maker <- function (ld, row, col) {
    # the distance between the snps
    dist <- abs(chrom$pos[row] - chrom$pos[col])
    # we scale the ld_floor up based on how close snps are to each other, we
    # do this quadratically because of the quadratic relationship between
    # marker distance and ld
    ld_scaled <- ld_ceiling - (ld_ceiling - ld_floor) * (1 - (0.8 / dist))
    # print("ld_scaled")
    if (ld_scaled > ld_ceiling) {
      ld_scaled <- ld_ceiling
    }
    # if ld is between ld_scaled and ld_ceiling we say the markers are linked
    # and place an edge between them in the adjacency matrix
    ifelse(ld >= ld_scaled && ld <= ld_ceiling, 1, 0)
  }
  # the proper row and column names for mapply so that the proper distance
  # between markers can be calculated
  col_names <- as.integer(colnames(sub_mat))
  # print("col names")
  # print(col_names)
  row_names <- as.integer(rownames(sub_mat))
  # turn these names into indices so that we can get the proper distance
  # between markers in adj_vec_maker
  row_index <- as.integer(as.vector(matrix(rep(row_names,
    nrow(sub_mat)), nrow(sub_mat), nrow(sub_mat))))
  # create a vector of the col name of each snp in sub_mat
  col_index <- as.integer(as.vector(matrix(rep(col_names,
    ncol(sub_mat)), ncol(sub_mat), ncol(sub_mat), byrow = TRUE)))
  # maplly returns a vector of the adjancency matrix with the values in the
  # order of concetenated columns
  adj_vec <- mapply(adj_vec_maker, sub_mat, row_index, col_index)
  # we build the adjacency matrix from the vector column first, dimnames are
  # those of the snp index value so that we can return the correct index values
  # below
  adj_mat <- matrix(adj_vec, nrow(sub_mat), ncol(sub_mat),
    dimnames = list(row_names, col_names))
  # finds the largest clique or cliques in the graph if two or more are of
  # equal length
  cliques <- graph_from_adjacency_matrix(adj_mat, mode = "undirected") %>%
    largest_cliques()
  # print("cliques")
  # print(cliques[[1]])
  # a vector to hold the clique members
  clique <- vector()
  # we intiailze the center of the clique as na which we can ignore if no
  # cliques over 3 members are found
  center <- NA
  # if the length of the largest clique(s) is over three members
  if (length(cliques) && length(cliques[[1]]) > 3) {
    # this pulls out the clique members with thier index values in ld_mat 
    # intact in sorted order as a vector
    clique <- cliques[[1]] %>% unlist() %>% as_ids() %>% as.integer() %>%
      sort()
    # this pulls out the clique mambers in index values of the sub_mat as a
    # vector in sorted order
    clique_pos <- cliques[[1]] %>% unlist() %>% as.integer() %>% sort()
    # finds the range of markers from the first_relative to last_relative in the clique, this
    # catches snps not in the clique but physically interlaced with it
    clique_pos_range <- clique_pos[1]:clique_pos[length(clique_pos)]
    # calcluates the average position of markers in the clique
    center <- mean(clique_pos_range)
    # subsets the submat to exclude all markers in the range of the clique
    sub_mat <- sub_mat[-clique_pos_range, -clique_pos_range]
  } else {
    # if no largest cliques > 3 are found we return an empty sub_mat
    sub_mat <- matrix()
  }
  # if no cliques were found clique is an empty vector and center is NA
  return (list(sub_mat = sub_mat, clique = clique, center = center))
}

find_first_last_relative <- function (densest_relative, window_size, pos,
  snp_index, genome_mean_dist) {
    start <- densest_relative
    end <- densest_relative
    max_dist <- (window_size + 1) * genome_mean_dist
    first_relative <- 1
    last_relative <- 1
    for (i in 1:(window_size %/% 2)) {
      dist <- pos[end] - pos[start]
      if (start - i >= 1 &&
        end + i <= length(snp_index) && dist < max_dist) {
          start <- start - i
          last_relative <- end + i
      } else if (start - i >= 1 && dist < max_dist) {
        start <- start - i
      } else if (end + i <= length(snp_index) && dist < max_dist) {
        end <- end + 1
      } else {
        first_relative <- start
        last_relative <- end
      }
    }
    return (list(first_relative = first_relative, 
      last_relative = last_relative))
  }

kept_finder <- function(chrom, genome_mean_dist, ld_mat, window_size,
  ld_floor, ld_ceiling, removed, remain) {
    print(3)
    # print(removed)
    snp_index <- 1:nrow(chrom)
    if (length(removed)) {
      chrom <- chrom[-removed, ]
      ld_mat <- ld_mat[-removed, -removed]
      snp_index <- snp_index[-removed]
    }
    # find the index of the first_relative marker with the max density in the
    # subset of remaining markers
    densest_relative <- which.max(snp_densities(chrom))
    # the first_relative and last_relative snps of the window centered on the
    # densest relative snp the max window size is 1/2 of window size upstream
    # and downstream of the densest snp as long as they are within window
    # size + 1 times the genome average marker distance apart
    first_last_relative <- find_first_last_relative(densest_relative,
      window_size, chrom$pos, snp_index, genome_mean_dist)
    first_relative <- first_last_relative$first_relative
    last_relative <- first_last_relative$last_relative
    # print("first_relative snp")
    # print(first_relative)
    # print("last_relative snp")
    # print(last_relative)
    # print("nrow ld_mat")
    # print(nrow(ld_mat))
    snp_range <- last_relative - first_relative + 1
    # print("snp_range")
    # print(snp_range)
    # print("first")
    # print(first_relative)
    # print("last")
    # print(last_relative)
    # the matrix of lds between markers in the window
    sub_mat <- matrix(ld_mat[first_relative:last_relative,
      first_relative:last_relative],
      snp_range, snp_range, dimnames = list(
        first_relative:last_relative, first_relative:last_relative))
    # print("nrow submat")
    # print(nrow(sub_mat))
    # print("sub mat")
    # print(sub_mat)
    # a temporary vector to contain the snps encompassed by the largest cliques
    # in the sub mat, including those in the clique
    removed_temp <- vector()
    # a vector to contain the physical centers of each clique
    centers <- vector()
    # find the first_relative largest clique in the sub matrix
    result <- largest_clique_finder(chrom, sub_mat, ld_floor, ld_ceiling)
    kept <- vector()
    # if result clique has a length greater than 0 process the result
    while (length(result$clique) && is.matrix(result$sub_mat) &&
      nrow(result$sub_mat)) {
        # print("reult sub mat")
        # print(result$sub_mat)
        # print("removed")
        # print(sort(removed))
        # print("clique members")
        # print(snp_index[result$clique])
        # add those snps in a clique retained from the window to the kept list
        kept <- c(kept, snp_index[result$clique])
        # those snps from the most downstream snp to the most upstream snp
        # encompassed by the range of the clique
        removed_actual <-
          snp_index[result$clique[1]:result$clique[length(result$clique)]]
        removed_temp <- c(removed_temp, removed_actual)
        # if a clique was returned the center will not be NA, we add this 
        # center to the vector of centers
        centers <- c(centers, result$center)
        # find the next largest clique in the window greater than 3
        result <- largest_clique_finder(chrom, result$sub_mat, ld_floor,
          ld_ceiling)
    }
    # if any snps were removed
    if (length(removed_temp)) {
      # print("removed_temp")
      # print(removed_temp)
      # indicate that these snps have been observed in the reamin vector
      remain[removed_temp] <- FALSE
      # add these snps to the removed list
      removed <- unique(c(removed, removed_temp))
    } else {
      # if there were no cliques found around the marker in its window indicate
      # that the markers has been observed and add it ro the removed vector
      # should i indicate that all snps in the window have been observed?
      densest_actual <- snp_index[densest_relative]
      remain[densest_actual] <- FALSE
      removed <- c(removed, densest_actual)
    }
    # print("removed after")
    # print(sort(removed))
    return (list(kept = kept, removed = removed, centers = centers,
      remain = remain))
}

ld_prune <- function (wheat, maf, window_size, ld_floor, ld_ceiling) {
  # should implement a method that finds the distance where ld equals
  # ld_ceiling on average and use that as the numerator in the scaling function
  print(1)
  # from the gds file pull out the information on the snps
  snp_data <- load_data_internal(wheat, maf)
  # calcualte the mean distance between snps on the same chromosome
  genome_mean_dist <- calc_genome_mean_dist(snp_data)
  # for each chromosome find snps to keep
  by(snp_data, snp_data$chrom, function (chrom) {
    # calcualte the ld matrix
    ld_mat <- abs(snpgdsLDMat(
        wheat, method = "composite", snp.id = chrom$id, slide = -1
    )$LD)
    # create a vector to contain the snps we are keeping
    kept <- vector()
    # create a vector to contain the snps we have removed from the data set
    # this includes kept snps and snps not in a clique but between the first_relative
    # and last_relative marker in a clique on the chromosome
    removed <- vector()
    # a vector to hold the physical middle positons of cliques, will be used
    # to calculate the average distance between clusters so that we can find
    # the amount of genome is covered by the cliques based on how much of the
    # genome is not covered by the range of 1/2 upstream and downstream of the
    # average distance from each clique center, will be used as the metric to
    # determine the qulaity of the pruning, is not sensitive to differential
    # ld between regions which may be problematic
    centers <- vector()
    # a vector containg the truth on whether a snp has been observed yet true
    # indicates the snp remains to be observed, false indicates the marker has
    # been found in a clique, between the start and end snp of a clique or does
    # not form a clique > n members with any nearby markers
    remain <- rep(TRUE, nrow(chrom))
    # find cliques until every marker on the chromosome has been observed
    result <- kept_finder(chrom, genome_mean_dist, ld_mat, window_size,
      ld_floor, ld_ceiling, removed, remain)
    # while any snps remain to be observed
    while (any(remain)) {
      # add the kept snps to the kept vector
      kept <- c(kept, result$kept)
      print("kept")
      print(kept)
      # add the removed snps to the reomoved vector
      removed <- result$removed
      # add the centers to the centers vector
      centers <- c(centers, result$centers)
      # update the remain vector to the returned one
      remain <- result$remain
      # pass in the data with the removed snps removed from the data set
      result <- kept_finder(chrom, genome_mean_dist, ld_mat, window_size,
        ld_floor, ld_ceiling, removed, remain)
    }
    # calculate the coverage statistic
    coverage <- centers %>%
      find_clique_regions() %>%
      calc_coverage()
    # return the kept snps and the coverage percentage of the chromosome
    return (list(ids = chrom$id[sort(kept)],
      coverage = coverage / chrom$pos[length(chrom$pos)]))
  })
}