phi.all <- function (wheat.amova) {
  phi.all <- vector()
  for (i in 1:length(wheat.amova)) {
    if (is.na(wheat.amova[[i]])) {
      phi.all[i] <- 0
    } else {
      phi.all[i] <- abs(wheat.amova[[i]]$statphi[1,1])
    }
  }
  return(phi.all)
}

sw_calc <- function (pos, value, ave = TRUE) {
  baseWidth <- mean(594102056, 689851870, 495453186, 780798557, 801256715, 651852609,
                    750843639, 830829764, 615552423, 744588157, 673617499, 509857067,
                    709773743, 713149757, 566080677, 618079260, 720988478, 473592718,
                    736706236, 750620385, 638686055)/200
  windowNear <- 0
  windowFar <- (baseWidth * 5) - 1
  
  i = 1
  meanValues <- vector()
  meanPoss <- vector()
  while (windowNear <= max(pos)) {
    index <- which(pos > windowNear & pos < windowFar)
    if (length(index) >= 10) {
      valueSet <- value[index]
      posSet <- pos[index]
      windowWidth <- 4
      if (ave) {
        for (j in 1:(length(posSet) - windowWidth)) {
          meanPoss[i] <- mean(posSet[j : (j + windowWidth)], na.rm = T)
          meanValues[i] <- mean(valueSet[j : (j + windowWidth)], na.rm = T)
          i = i + 1
        }
      } else {
        for (j in 1:(length(posSet) - windowWidth)) {
          meanPoss[i] <- sd(posSet[j : (j + windowWidth)], na.rm = T)
          meanValues[i] <- sd(valueSet[j : (j + windowWidth)], na.rm = T)
          i = i + 1
        }
      }
    } else {
      meanValues[i] <- 0
      meanPoss[i] <- mean(c(windowNear, windowFar))
      i = i + 1
    }
    windowNear = windowNear + baseWidth
    windowFar = windowFar + baseWidth
  }
  
  for (j in 1:length(meanValues)) {
    if (is.na(meanValues[j])) {
      meanValues[j] <- 0
    }
  }
  
  res <- cbind(meanPoss, meanValues)
  res <- res[order(res[,1]),]
  return(res)
}

EH <- function (genotypes) {
  apply(genotypes, 1, function (row) {
    2*((sum(row == 0)/sum(row == 0 | 2)) * (sum(row == 2)/sum(row == 0 | 2)))
  })
}