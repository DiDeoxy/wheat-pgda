library(SNPRelate)
library(igraph)
library(tidyverse)

parse_ld <- function (adj_mat, kept, start_snp, end_snp, pruned) {
  # print(7)
  # create the cliques for the window determined by start_snp to end_snp
  # (the max position in the window determined by end_snps[start_snp])
  cliques_frame <- cliques_frame_builder(start_snp, end_snp, adj_mat)
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
  # if there are no cliques and no end_snps in the window have been kept
  # or are in the pruned vector add end_snp to the end of the list
  } else if (! any(start_snp:end_snp %in% kept) &&
    ! any(start_snp:end_snp %in% pruned)) {
    return (c(kept, end_snp))
  }
  # return the kept list
  kept
}

cliques_frame_maker <- function (cliques) {
  # print(6)
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

cliques_frame_builder <- function (start_snp, end_snp, adj_mat) {
  # print(5)
  # subset the adjacency matrix to just those from start_snp to end_snp
  # (max position in num_snps)
  sub_mat <- matrix(
    adj_mat[start_snp:end_snp, start_snp:end_snp],
      (end_snp - start_snp + 1), (end_snp - start_snp + 1),
      dimnames = list(start_snp:end_snp, start_snp:end_snp))
  # find the maximal cliques with a min size determined by the number of
  # snps between start_snp and end_snp and turn these into a
  # data frame for easier parsing and return the cliques frame
  graph_from_adjacency_matrix(sub_mat, mode = "undirected") %>%
    max_cliques(min = ((end_snp - start_snp) %/% 3)) %>%
    cliques_frame_maker()
}

adj_mat_perfect_rep_maker <- function (ld_mat, end_snps) {
  # print(4)
  # create an adjacency matrix where all snpsin perfect ld ar connected
  adj_mat_perfect <- ifelse(ld_mat > 0.99 & ld_mat < 1.01, 1, 0)
  # since all snps compared to themself are in perfect ld set these values
  # to zero to remove the loop
  for (start_snp in 1:nrow(adj_mat_perfect)) {
    adj_mat_perfect[start_snp, start_snp] <- 0
  }
  # create an empty adjacenecy matrix
  adj_mat_perfect_rep <- matrix(rep(0, length(adj_mat_perfect)),
    nrow(adj_mat_perfect), ncol(adj_mat_perfect))
  # a vector to contain all snps not sample from the cliques
  pruned <- vector()
  all_cliques <- vector()
  for (start_snp in 1:nrow(adj_mat_perfect)) {
    # creat a data from of all the cliques
    cliques_frame <- cliques_frame_builder(start_snp,
      end_snps[start_snp], adj_mat_perfect)
    if (nrow(cliques_frame)) {
      for (start_snp in 1:nrow(cliques_frame)) {
        # find all end_snps in the clique that aren't na
        clique <- cliques_frame[start_snp,
          ! is.na(cliques_frame[start_snp, ])]
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
          for (snp in 1:ncol(pairs)) {
            combo <- unlist(pairs[, snp])
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

adj_mat_range_maker <- function (pos, ld_mat, ld_floor, ld_ceiling,
  end_snps) {
    # print(3)
    # create matrix of the same size as ld_mat that has only zeroes
    adj_mat_range <- matrix(rep(0, length(ld_mat)), nrow = nrow(ld_mat),
      ncol = ncol(ld_mat))
    # for each row in ld_mat
    for (start_snp in 1:(nrow(ld_mat) - 1)) {
      # find the distance between snp start_snp and snp end_snp for
      # the row
      max_dist <- pos[end_snps[start_snp]] - pos[start_snp]
      i_row <- vector()
      # for each snp from start_snp + 1 to start_snp + n
      for (snp in (start_snp + 1):end_snps[start_snp]) {
        # find the distance between snp snp of and snp start_snp
        dist <- pos[snp] - pos[start_snp]
        # find the ld between these snps
        ld <- ld_mat[start_snp, snp]
        # scale the ld_floor so that closer end_snps have to have higher ld
        # to each other to be included
        ld_scaled <- ld_ceiling - ((ld_ceiling - ld_floor) *
          (1 - (dist / max_dist) ^ 2))
        # if ld between the snps is between ld_scaled and ld_ceiling place
        # a 1 in i_row indicating the snps are connected
        i_row <- c(i_row, ifelse(ld >= ld_scaled && ld <= ld_ceiling, 1, 0))
      }
      # add the snps in i_row to the adj cency matrix
      adj_mat_range[start_snp,
        (start_snp + 1):end_snps[start_snp]] <- i_row
      adj_mat_range[(start_snp + 1):end_snps[start_snp],
        start_snp] <- i_row
    }
    # return the adjacency matrix
    adj_mat_range
  }

find_snps <- function (pos, min_snps, max_snps, num_rows, genome_mean_dist) {
  # print(2)
  end_snps <- vector()
  # for every row except the last one
  for (start_snp in 1:(num_rows - 1)) {
    # find the middle value between max_snps and min_snps
    mean_snps <- mean(c(max_snps, min_snps))
    # the potential max number of snps in num_snps is max_snps
    # unless we are getting close to the end of a chromosome
    max_snps <- ifelse(start_snp + max_snps <= num_rows,
      max_snps, num_rows - start_snp)
    # find the mean distance between snps from start_snp to 
    # start_snp + max_snps
    local_mean_dist <- mean(diff(pos[start_snp:(start_snp + max_snps)]))
    # find the log10 of the ratio of the genome wide mean snp distance
    # over the local mean snp distance
    log_ratio_genome_local <- log10(genome_mean_dist / local_mean_dist)
    # num_snps is the number of snps in the window, it is the
    # difference between the max and mean number of snps times the log
    # ratio plus mean_snps, this means that denser regions will
    # have more snps in their windows than less dense regions
    num_snps <- round(mean_snps + (max_snps - mean_snps) *
      log_ratio_genome_local)
    # join start_snp + num_snps to snps if num_snps is between
    # min_snps and max_snps and start_snp + num_snps is less than
    # num_rows, else depending on if num_snps is less than min_snps or
    # greater than max_snps add start_snp + those, else return num_rows
    if (num_snps <= max_snps && num_snps >= min_snps &&
      start_snp + num_snps <= num_rows) {
        end_snps <- c(end_snps, start_snp + num_snps)
    } else if (num_snps >= max_snps &&
      start_snp + max_snps <= num_rows) {
        end_snps <- c(end_snps, start_snp + max_snps)
    } else if (num_snps <= min_snps &&
      start_snp + min_snps <= num_rows) {
        end_snps <- c(end_snps, start_snp + min_snps)
    } else {
        end_snps <- c(end_snps, num_rows)
    }
  }
  # return the end_snps plus an extra num_rows for the final row
  c(end_snps, num_rows)
}

adj_mat_maker <- function (pos, min_snps, max_snps, ld_mat, ld_floor,
  ld_ceiling, genome_mean_dist) {
    # print(1)
    # for each start_snp find end_snp which is the downstream snp
    # fulfilling the search criteria of find_snps
    end_snps <- find_snps(pos, min_snps, max_snps, nrow(ld_mat),
      genome_mean_dist)
    # in the case where ld_ceiling is below the max possible ld
    if (ld_ceiling < 1) {
      # find adjacency between snps with ld between ld_floor and ld_ceiling
      adj_mat_range <- adj_mat_range_maker(pos, ld_mat, ld_floor, ld_ceiling,
        end_snps)
      # find the adjacency matrix of samples of cliques of snps in perfect
      # ld to each other and the set of snps excluded from these samples
      temp <- adj_mat_perfect_rep_maker(ld_mat, end_snps)
      # add the two adjacency matrices together to create the final set
      adj_mat_both <- adj_mat_range + temp$adj_mat_perfect_rep
      # all snps that are in both matrices will have a value of two, we
      # change this to 1, this should be rare
      adj_mat_both[adj_mat_both == 2] <- 1
      # return the adjacency matrix, the pruned set of snps, and
      # end_snps
      return (list(adj_mat = adj_mat_both, pruned = temp$pruned,
        end_snps = end_snps))
    } else {
      # in the case where the ld_ceiling is above equal to or greater than 1
      # we can skip making adj_mat_perfect_rep
      adj_mat_range <- adj_mat_range_maker(pos, ld_mat, ld_floor, ld_ceiling,
        end_snps)
      # there are no pruned snps so we return an empty vector
      return (list(adj_mat = adj_mat_range, pruned = vector(),
        end_snps = end_snps))
    }
  }

ld_prune <- function (wheat, maf, min_snps, max_snps, ld_floor,
  ld_ceiling) {
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
    # ld prune the end_snps on each chromosome
    by(snp_data, snp_data$chrom, function (chrom) {
      # construct a matrix of the ld between all pairs of snps on the chrom
      ld_mat <- abs(snpgdsLDMat(
          wheat, method = "composite", snp.id = chrom$id, slide = -1
      )$LD)
      # create an adjancency matrix using the ld_matrix
      temp <- adj_mat_maker(chrom$pos, min_snps, max_snps, ld_mat,
        ld_floor, ld_ceiling, genome_mean_dist)
      # create an empty list to hold the snps we keep
      kept <- vector()
      # for each snp in chrom
      for (start_snp in 1:nrow(chrom)) {
        # find the snps in cliques not previously kept or in the pruned list
        kept <- parse_ld(temp$adj_mat, kept, start_snp, 
          temp$end_snps[start_snp], temp$pruned)
      }
      # return the ids of the kept snps
      chrom$id[sort(kept)]
    })
  }