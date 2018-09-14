phi_markers <- function(amova) {
  phis <- vector()
  for (i in 1:length(amova)) {
    if (length(amova[[i]]) == 0 || is.na(amova[[i]])) {
      phis[i] <- 0
    } else {
      phis[i] <- abs(amova[[i]]$statphi[1, 1])
    }
  }
  return(phis)
}

eh <- function(genotypes) {
  apply(genotypes, 1, function(row) {
    2 * ((sum(row == 0) / sum(row == 0 | 2)) *
      (sum(row == 2) / sum(row == 0 | 2)))
  })
}