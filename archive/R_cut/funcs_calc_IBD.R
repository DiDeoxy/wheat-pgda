windowMaker <- function (pos, genotypes.imputed, chroms) {
  windowsPos <- by(cbind(pos, t(genotypes.imputed[,-1])), chroms, function (chrom) {
    pos <- chrom[,1]
    chrom <- t(as.matrix(chrom[,-1]))
    # windowsPos <- list
    windows <- apply(chrom, 1, function (indiv) {
      windowsSnps <- list()
      windowsPos <- list()
      i <- 1
      sapply(seq(1, length(indiv), 15), function (windowStart) {
        if (windowStart + 19 <= length(indiv)) {
          windowEnd <- windowStart + 19
          windowsSnps[[i]] <<- indiv[windowStart:windowEnd]
          windowsPos[[i]] <<- c(pos[windowStart], pos[windowEnd])
        } else {
          if (length(indiv) >= 20) {
            windowStart <- length(indiv) - 19
          }
          windowEnd <- length(indiv)
          windowsSnps[[i]] <<- indiv[windowStart:windowEnd]
          windowsPos[[i]] <<- c(pos[windowStart], pos[windowEnd])
        }
        i <<- i + 1
      })
      list(windowsSnps, windowsPos)
    })
  })
  return(windowsPos)
}

hapLengthsTotalFinder <- function (windowsPos) {
  chromsDistances <- lapply(windowsPos, function (chrom) {
    chromDistances <- matrix(nrow = 358, ncol = 358, data = 0)
    for (i in 1:length(chrom)) {
      indivISegs <- chrom[[i]][[1]]
      indivIDists <- chrom[[i]][[2]]
      for (j in i+1:length(chrom)) if (j < length(chrom)) {
        indivJSegs <- chrom[[j]][[1]]
        start <- indivIDists[[1]][1]
        end <- indivIDists[[1]][2]
        lastIdentity <- sum(indivISegs[[1]] == indivJSegs[[1]], na.rm = T) / length(indivISegs[[1]])
        for (segment in 1:length(indivISegs)) {
          segIdentity <- sum(indivISegs[[segment]] == indivJSegs[[segment]], na.rm = T) / length(indivISegs[[segment]])
          if (segIdentity >= 0.90 & lastIdentity >= 0.90) {
            end <- indivIDists[[segment]][2]
          } else if (segIdentity <= 0.90 & lastIdentity >= 0.90) {
            chromDistances[i,j] <- chromDistances[i,j] + (end - start)
            if (segment + 1 <= length(indivISegs)) {
              start <- indivIDists[[segment + 1]][1]
              end <- indivIDists[[segment + 1]][2]
            }
          } else {
            if (segment + 1 <= length(indivISegs)) {
              start <- indivIDists[[segment + 1]][1]
              end <- indivIDists[[segment + 1]][2]
            }
          }
          lastIdentity <- segIdentity
        }
      }
    }
    chromDistances
  })
  
  return(totalDistances <- Reduce('+', chromsDistances))
}