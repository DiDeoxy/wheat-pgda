allele.richness2 <- function (pop, coding) {
  n <- ncol(pop)
  markers <- rowTables(pop, coding)
  lapply(0:(n - 1), function (k) {
    lapply(1:nrow(markers), function(marker) {
      (1 - lapply(1:length(markers[marker, ]), function (allele) {
        lapply(0:k, function (u) {
          val <- (n - markers[[marker, allele]] - u) / (n - u)
          if (val > 0) {
            val
          } else {
            0
          }
        }) %>% unlist() %>% prod()
      }) %>% unlist()) %>% sum()
    }) %>% unlist() %>% mean()
  }) %>% unlist()
}