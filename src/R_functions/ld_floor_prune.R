ld_floor_prune <- function(chrom, ld_mat, ld_floor, max_window) {
  # pick a random starting point on each chromosome
  start <- sample(1:nrow(chrom), 1)
  # start the sliding window at this point
  in_window <- start
  kept <- vector()
  for (i in (start + 1):nrow(chrom)) {
    temp <- slide_window(
      chrom, kept, in_window, i, max_window, ld_mat, ld_floor
    )
    if (length(temp) == 1 && ! temp) {
      next
    } else if (length(temp == 1)) {
      in_window <- temp[[1]]
    } else {
      kept <- temp[1]
      in_window <- temp[2]
    }
  }
  in_window <- start
  for (i in (start - 1):1) {
    temp <- slide_window(
      chrom, kept, in_window, i, max_window, ld_mat, ld_floor
    )
    if (length(temp) == 1 && ! temp) {
      next
    } else if (length(temp == 1)) {
      in_window <- temp[[1]]
    } else {
      kept <- temp[1]
      in_window <- temp[2]
    }
  }
  return (kept)
}

slide_window <- function (chrom, kept, in_window, i, max_window, ld_mat,
  ld_floor) {
    # prune the window so that only markers less than max_window distant from
    # the ith marker or only one marker remains in the window, pruned markers
    # are added to the kept list
    temp <- prune_window(chrom, kept, in_window, i, max_window)
    if (length(temp) == 1) {
      in_window <- temp
    } else {
      kept <- temp[1]
      in_window <- temp[2]
    }
    # if window has only one marker in it and the position of the next marker
    # is further away than the max_window add the current marker to the kept
    # list and make the next marker the only marker in the window, skip to the
    # next iteration of the loop
    if (length(in_window) == 1)
      if (abs(chrom$pos[i] - chrom$pos[in_window]) > max_window) {
        kept <- in_window
        in_window <- i
        return(c(kept, in_window))
      # there is a possibility if we have only one marker in the window and the
      # next marker is within max_window that the marker in the window is not
      # from that locus. To overcome this we test all combinations of markers
      # within max_window until we find a pair in ld above ld_floor and make
      # these the markers in the sliding window dropping the one that was
      # already there and moving onto the next iteration of the loop
      } else if (ld_mat[in_window, i] < ld_floor) {
        j <- 1
        while (abs(chrom$pos[i + j] - chrom$pos[in_window]) < max_window) {
          j <- j + 1
        }
        combos <- combn(in_window:(i + j))
        for (k in ncol(combos)) {
          if (ld_mat[combos[1, k], combo[2, k]] > ld_floor) {
            return(list(combos[, k]))
          }
        }
        # if none of the combos above found markers in ld above the floor we
        # drop all markers from the locus and jump ahead to the (i + j)th
        # marker
        if (length(in_window) == 1) {
          in_window <- i + j
          return(0)
        }
      }
    # since we potentially have marerks in the window greater than i due to
    # above step we check the last marker in the window and skip ahead until
    # the ith marker is greater
    if (in_window[length(in_window)] >= i) {
      return(0)
    }
    # the ith marker is guranteed to be within max_window to the first marker
    # in the window, if it is above the ld_floor to any marker in the window
    # add it to the window
    in_window <- c(in_window,
      check_floor(in_window, i, ld_mat, ld_floor)
    )
    return(c(kept, in_window))
  }

prune_window <- function (chrom, kept, in_window, i, max_window) {
  # if the window has more than one marker in it recursively remove markers and
  # add them to the kept list until one marker is left in the window or the
  # ith marker is less than max_window distant from the first marker in the
  # window
  if (length(in_window) > 1 &&
    abs(chrom$pos[i] - chrom$pos[in_window[1]] > max_window)) {
      return(
        prune_window(
          chrom, c(kept, in_window[1]), in_window[2:length(in_window)], i,
          max_window
        )
      )
  } else {
    return(c(kept, in_window))
  }
}

check_floor <- function(in_window, i, ld_mat, ld_floor) {
  # if the ith marker has greater ld than ld_floor to the jth marker in the
  # sliding window add the ith marker to the sliding window
  for (j in in_window) {
    if (ld_mat[j, i] > ld_floor) {
        return (i)
    }
  }
  return (NULL)
}
