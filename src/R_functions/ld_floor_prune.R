ld_floor_prune <- function(chrom, ld_mat, ld_floor, max_window) {
  # make sure ld_mat has only distances from zero in it
  ld_mat <- abs(ld_mat)
  # pick a random starting point on each chromosome
  start <- sample(1:nrow(chrom), 1)
  # initialize the kept list of markers with NA, to be removed later
  kept <- NA
  # window heading upstream
  in_window <- start
  for (i in ((start + 1):nrow(chrom))) {
    temp <- slide_window(
      1, chrom, kept, in_window, i, max_window, ld_mat, ld_floor
    )
    kept <- temp[[1]]
    in_window <- temp[[2]]
  }
  # window heading downstream
  in_window <- start
  for (i in (start - 1):1) {
    temp <- slide_window(
      -1, chrom, kept, in_window, i, max_window, ld_mat, ld_floor
    )
    kept <- temp[[1]]
    in_window <- temp[[2]]
  }
  return (na.omit(kept))
}

slide_window <- function (dir, chrom, kept, in_window, i, max_window, ld_mat,
  ld_floor) {
    # we potentially have markers in the window greater than i due to a
    # proceeding step. we therefor check if the ith marker is before the the
    # last marker in the window, if it is we return the kept and in_window 
    # vectors as is
    if (dir > 0 && i < in_window[length(in_window)]) {
      # print("1a")
      return(list(kept, in_window))
    } else if (dir < 0 && i > in_window[length(in_window)]) {
      # print("1b")
      return(list(kept, in_window))
    }
    # prune the window so that only markers less than max_window distant from
    # the ith marker or only one marker remains in the window, pruned markers
    # are added to the kept list
    temp <- prune_window(chrom, kept, in_window, i, max_window)
    kept <- temp[[1]]
    in_window <- temp[[2]]
    # if window has only one marker in it and the position of the ith marker
    # is further away than max_window retunr the concatenation of the kept
    # and in_window list as the new kept list and i as the start of the window
    if (length(in_window) == 1) {
      # if (abs(chrom$pos[i] - chrom$pos[in_window]) < max_window) {
      #   print(ld_mat[in_window, i])
      # }
      if (abs(chrom$pos[i] - chrom$pos[in_window]) > max_window) {
        # print(2)
        return(list(c(kept, in_window), i))
      # there is a possibility if we have only one marker in the window and the
      # ith marker is within max_window that the marker in the window is not
      # from that locus. To overcome this we test all combinations of markers
      # within max_window until we find the earliest pair in ld above ld_floor
      # we then return the kep list as is an make the sliding window the 
      # markers in ld  floor
      } else if (is.nan(ld_mat[in_window, i]) ||
        ld_mat[in_window, i] < ld_floor) {
        j <- 1
        while ((i + j) < nrow(chrom) && 
          abs(chrom$pos[i + j] - chrom$pos[in_window]) < max_window) {
            j <- j + 1
        }
        combos <- combn(in_window:(i + j), 2)
        for (k in 1:ncol(combos)) {
          if (! is.nan(ld_mat[combos[1, k], combos[2, k]]) && 
            ld_mat[combos[1, k], combos[2, k]] > ld_floor) {
            # print(3)
            return(list(kept, combos[, k]))
          }
        }
        # if none of the combos above found markers in ld above the floor we
        # return the kept list as is and make the (i + j)th  marker the start 
        # of the window
        # print(4)
        return(list(kept, i + j))
        # i think this would be worse:
        # set.seed(1000)
        # return(c(c(kept, sample(i:j, 1)), i + j)) 
      # since the distance to ith marker from the marker in the window is less 
      # than max_window and the ld between the two markers is above ld_floor
      # return the kept list as is and add the ith marker to the window
      } else {
        # print(5)
        return(list(kept, c(in_window, i)))
      }
    # if length(in_window) > 1, the ith marker will always be less than 
    # max_window distant because of the pruning step, thus we test every marker
    # j in the window against the ith marker, if any have ld above ld_floor
    # the ith marker is added to the window
    } else {
      for (j in in_window) {
        if (! is.nan(ld_mat[j, i]) && ld_mat[j, i] > ld_floor) {
          # print(6)
          return(list(kept, c(in_window, i)))
        }
      }
      # if the ith marker is below ld_floor to all markers in the window return
      # kept and in_window as is
      return(list(kept, in_window))
    }
    # print("no")
  }

prune_window <- function (chrom, kept, in_window, i, max_window) {
  # if the window has more than one marker in it recursively remove markers and
  # add them to the kept list until one marker is left in the window or the
  # ith marker is less than max_window distant from the first marker in the
  # window
  if (length(in_window) > 1 &&
    abs(chrom$pos[i] - chrom$pos[in_window[1]]) > max_window) {
      return(
        prune_window(
          chrom, c(kept, in_window[1]), in_window[2:length(in_window)], i,
          max_window
        )
      )
  } else {
    # print("yes")
    return(list(kept, in_window))
  }
}
