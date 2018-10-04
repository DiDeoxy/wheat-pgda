library(SNPRelate)
library(igraph)
library(tidyverse)

parse_ld <- function (chrom, adj_mat, kept, i, j, pruned) {
  # create the cliques for the max_snps determined by i to j (the max position
  # in the max_snps determined by windows[i])
  cliques_frame <- cliques_frame_builder(i, j, chrom, adj_mat)
  if (nrow(cliques_frame)) {
    for (clique in cliques_frame) {
      # remove nas in the clique row
      clique <- unlist(na.omit(clique))
      for (snp in clique) {
        # for each snp in the clique check that it is not already being kept
        # and that it wasn't pruned in the construction of adj_mat_perfect_rep
        # if not add it to the kept list
        if (! snp %in% kept && ! snp %in% pruned) {
          kept <- c(kept, snp)
        }
      }
    }
  # if there are no cliques and no snps in the max_snps have been kept or are
  # in the pruned vector add j to the end of the list
  } else if (! any(i:j %in% kept) && ! any(i:j %in% pruned)) {
    return (c(kept, j))
  }
  # return the kept list
  kept
}

cliques_frame_maker <- function (cliques) {
  # take cliques and turn them into a data frame, shorter cliques have their
  # missing values filled in with NA
  cliques_frame <- tibble()
  if (length(cliques)) {
    lengths <- lapply(cliques, length)
    max_length <- lengths[[which.max(lengths)]]
    for (clique in cliques) {
      clique <- clique %>%
        as_ids() %>%
        as.integer() %>%
        sort() %>%
        unlist()
      length(clique) <- max_length
      cliques_frame <- rbind(cliques_frame, clique)
    }
  }
  # return the cliques frame
  cliques_frame
}

cliques_frame_builder <- function (i, j, chrom, adj_mat) {
  # subset the adjacency matrix to just those from i to j (max position in
  # max_snps)
  sub_mat <- matrix(adj_mat[i:j, i:j], (j - i + 1), (j - i + 1),
    dimnames = list(i:j, i:j))
  # find the maximal cliques with a min size determined by the number of
  # snps between i and j and turn these into a data frame for easier\
  # parsing and return the cliques frame
  graph_from_adjacency_matrix(sub_mat, mode = "undirected") %>%
    max_cliques(min = ((j - i) %/% 2)) %>%
    cliques_frame_maker()
}

adj_mat_perfect_rep_maker <- function (chrom, ld_mat, windows) {
  # create an adjacency matrix where all snps in perfect ld are connected
  adj_mat_perfect <- ifelse(ld_mat > 0.99 & ld_mat < 1.01, 1, 0)
  # since all snps compared to themself are in perfect ld set these values 
  # to zero to remove the loop
  for (i in 1:nrow(adj_mat_perfect)) {
    adj_mat_perfect[i, i] <- 0
  }
  # create an empty adjacenecy matrix
  adj_mat_perfect_rep <- matrix(rep(0, length(adj_mat_perfect)),
    nrow(adj_mat_perfect), ncol(adj_mat_perfect))
  # a vector to contain all snps not sample from the cliques
  pruned <- vector()
  all_cliques <- vector()
  for (i in 1:nrow(adj_mat_perfect)) {
    # creat a data from of all the cliques
    cliques_frame <- cliques_frame_builder(i, windows[i], chrom,
      adj_mat_perfect)
    if (nrow(cliques_frame)) {
      for (i in 1:nrow(cliques_frame)) {
        # find all snps in the clique that aren't na
        clique <- cliques_frame[i, ! is.na(cliques_frame[i, ])]
        # if there are atleast 3 makers in the clique and no snps are in 
        # the pruned set or from previous cliques
        if (length(clique) >= 3 && ! any(clique %in% pruned) &&
          ! any(clique %in% all_cliques)) {
          # sample 3 snps from the clique
          set.seed(1000)
          clique_sample_index <- sample(1:length(clique), 3)
          # pull out the sample snps
          clique_sample <- clique[clique_sample_index]
          # add the snps not in the sample to the pruned set, this will be
          # returned so that snps in it are not included in the set of
          # kept snps from the chromosome
          pruned <- c(pruned, clique_sample[-clique_sample_index])
          # store all snps from the clique, this will be used internally to
          # prevent cliques with snps from another clique from being
          # included in the adjacency matrix
          all_cliques <- c(all_cliques, clique)
          # find all pairs of the three snps
          pairs <- combn(clique_sample, 2)
          # add an edge between each pair of snps
          for (k in 1:ncol(pairs)) {
            combo <- unlist(pairs[, k])
            adj_mat_perfect_rep[combo[1], combo[2]] <- 1
            adj_mat_perfect_rep[combo[2], combo[1]] <- 1
          }
        }
      }
    }
  }
  # return the adjacency matrix and the pruned set of snps
  return (list(adj_mat_perfect_rep = adj_mat_perfect_rep, pruned = pruned))
}

adj_mat_range_maker <- function (chrom, ld_mat, ld_floor, ld_ceiling, 
  windows) {
    # create matrix of the same size as ld_mat that has only zeroes
    adj_mat_range <- matrix(rep(0, length(ld_mat)), nrow = nrow(ld_mat), 
      ncol = ncol(ld_mat))
    # for each row in ld_mat
    for (i in 1:(nrow(ld_mat) - 1)) {
      # find the distance between snp i and snp j for the row
      max_dist <- chrom$pos[windows[i]] - chrom$pos[i]
      i_row <- vector()
      # for each snp from i + 1 to i + n
      for (k in (i + 1):windows[i]) {
        # find the distance between snp k of and snp i
        dist <- chrom$pos[k] - chrom$pos[i]
        # find the ld between these snps
        ld <- ld_mat[i, k]
        # scale the ld_floor so that closer snps have to have higher ld
        # to each other to be included
        ld_scaled <- ld_ceiling - ((ld_ceiling - ld_floor) *
          (1 - (dist / max_dist) ^ 2))
        # if ld between the snp is between ld_scaled and ld_ceiling place
        # a 1 in i_row indicating the snps are connected
        i_row <- c(i_row, ifelse(ld >= ld_scaled && ld <= ld_ceiling, 1, 0))
      }
      # add the snps in i_row to the adj cency matrix
      adj_mat_range[i, (i + 1):windows[i]] <- i_row
      adj_mat_range[(i + 1):windows[i], i] <- i_row
    }
    # return the adjacency matrix
    adj_mat_range
  }

find_windows <- function (chrom, max_snps, num_rows, genome_mean_dist) {
  windows <- vector()
  for (i in 1:num_rows) {
    # max_snps is max_snps if there is enough room else its the difference between
    # the last snp and n
    i_window <- ifelse(i + max_snps < num_rows, max_snps, num_rows - i)
    windows <- c(windows, i + i_window)
  }
  # return the windows
  windows
}

adj_mat_maker <- function (chrom, max_snps, ld_mat, ld_floor, ld_ceiling,
  genome_mean_dist) {
    # for each i find j which is the downstream snp fulfilling the search
    # criteria of find_windows
    windows <- find_windows(chrom, max_snps, nrow(ld_mat),
      genome_mean_dist)
    # in the case where ld_ceiling is below the max possible ld
    if (ld_ceiling < 1) {
      # find adjacency between snps with ld between ld_floor and ld_ceiling
      adj_mat_range <- adj_mat_range_maker(chrom, ld_mat, ld_floor, ld_ceiling,
        windows)
      # find the adjacency matrix of samples of cliques of snps in perfect
      # ld to each other and the set of snps excluded from these samples
      temp <- adj_mat_perfect_rep_maker(chrom, ld_mat, windows)
      # add the two adjacency matrices together to create the final set
      adj_mat_both <- adj_mat_range + temp$adj_mat_perfect_rep
      # all snps that are in both matrices will have a value of two, we
      # change this to 1, this should be rare
      adj_mat_both[adj_mat_both == 2] <- 1
      # return the adjacency matrix, the pruned set of snps, and
      # windows
      return (list(adj_mat = adj_mat_both, pruned = temp$pruned,
        windows = windows))
    } else {
      # in the case where the ld_ceiling is above equal to or greater than 1
      # we can skip making adj_mat_perfect_rep
      adj_mat_range <- adj_mat_range_maker(chrom, ld_mat, ld_floor, ld_ceiling,
        windows)
      # there are no pruned snps so we return an empty vector
      return (list(adj_mat = adj_mat_range, pruned = vector(),
        windows = windows))
    }
  }

ld_prune <- function (wheat, maf, max_snps, ld_floor, ld_ceiling) {
  # extract the snp metadata for the selected snp
  snp_id <- as.character(
    read.gdsn(index.gdsn(wheat, "snp.id")))
  snp_chrom <- as.integer(
    read.gdsn(index.gdsn(wheat, "snp.chromosome")))
  snp_pos <- as.integer(
    read.gdsn(index.gdsn(wheat, "snp.position")))
  # find the subset of snps with maf above the given threshold
  snp_maf <- snpgdsSelectSNP(wheat, autosome.only = F, maf = maf)
  snp_maf_indices <- match(unlist(snp_maf), snp_id)
  # create a tibble of the snp metadata
  snp_data <- tibble(
    id = snp_id[snp_maf_indices], chrom = snp_chrom[snp_maf_indices],
    pos = snp_pos[snp_maf_indices], stringsAsFactors = FALSE)
  # find the distances between neighbouring snps on each chromsomes
  dists <- by(snp_data, snp_data$chrom, function (chrom) {
    diff(chrom$pos)
  })
  # calc the average distance between snps genome wide
  genome_mean_dist <- mean(unlist(dists))
  # ld prune the snps on each chromosome
  by(snp_data, snp_data$chrom, function (chrom) {
    # construct a matrix of the ld between all pairs of snps on the chrom
    ld_mat <- abs(snpgdsLDMat(
        wheat, method = "composite", snp.id = chrom$id, slide = -1
    )$LD)
    # create an adjancency matrix using the ld_matrix
    temp <- adj_mat_maker(chrom, max_snps, ld_mat, ld_floor, ld_ceiling,
      genome_mean_dist)
    # create an empty list to hold the snps we keep
    kept <- vector()
    # for each snp in chrom
    for (i in 1:nrow(chrom)) {
      # k is greater than 0 decrment it by 1 and moce on to the next
      # iteration of the loop
      # find the snps in cliques not previously kept or in the pruned list
      kept <- parse_ld(chrom, temp$adj_mat, kept, i, temp$windows[i],
        temp$pruned)
    }
    # return the ids of the kept snps
    chrom$id[sort(kept)]
  })
}