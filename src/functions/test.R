calc_group_extreme_freqs <- function(wheat_data, extremes, prune = FALSE) {
  by(wheat_data$snp, wheat_data$snp$chrom,
    function(snp_data) {
      # initialize an empyt tibble with columns for the important data
      ret <- tibble(
        chrom = integer(), group = character(), pos_mb = double(),
        freq = double(), extreme_D = double(), id = character()
      )
      # 
      for (group in names(extremes)) {
        for (i in 1:length(snp_data$id)) {
        num_nearby <- which(
          snp_data$pos_mb >= snp_data$pos_mb[i] - 5 &
          snp_data$pos_mb <= snp_data$pos_mb[i] + 5
        ) %>% length()
        freq <- integer()
        if (num_nearby >= 5 && i >= 3 && i <= (length(snp_data$id) - 2)) {
          freq <- sum(
            snp_data[[group]][(i - 2):(i + 2)] > extremes[[group]]
          ) / 5
        } else if (
          (i < 3 || i > (length(snp_data$id) - 2)) &&
          snp_data[[group]][i] > extremes[[group]]
        ) {
          freq <- 0.1
        } else {
          freq <- NA
        }
        if (! is.na(freq) && freq == 0) {
          freq <- NA
        }
          ret <- ret %>%
            add_row(
              chrom = snp_data$chrom[i], group = group,
              pos_mb = snp_data$pos_mb[i], freq = freq,
              extreme_D = ifelse(
                snp_data[[group]][i] > extremes[[group]], snp_data[[group]][i],
                NA
              ),
              id = snp_data$id[i]
            )
        }
      }
      # set all markers that don't have the highest frequency of extreme nearby
      # markers to zero for each contiguous region of extreme markers
      if (prune) {
        temp <- vector()
        for (i in 1:nrow(ret)) {
          if (! is.na(ret[i, ]$freq)) {
            temp <- c(temp, i)
          } else if (length(temp) > 0) {
            highest <- which(ret[temp, "freq"][[1]] == max(ret[temp, "freq"]))
            ret[temp[-highest], "freq"] <- NA
            temp <- vector()
          }
        }
      }
      print(ret)
      ret
    }
  ) %>% do.call(rbind, .)
}
