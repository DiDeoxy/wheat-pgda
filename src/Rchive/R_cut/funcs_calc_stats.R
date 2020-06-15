phi_markers <- function (amova) {
  phis <- vector()
  for (i in 1:length(amova)) {
    if (length(amova[[i]]) == 0 || is.na(amova[[i]])) {
      phis[i] <- 0
    } else {
      phis[i] <- abs(amova[[i]]$statphi[1,1])
    }
  }
  return(phis)
}

# sw_calc <- function (pos, value, ave = TRUE) {
#   base_width <- mean(594102056, 689851870, 495453186, 780798557, 801256715,
#                      651852609, 750843639, 830829764, 615552423, 744588157,
#                     673617499, 509857067, 709773743, 713149757, 566080677,
#                     618079260, 720988478, 473592718, 736706236, 750620385,
#                     638686055) / 200
#   window_near <- 0
#   window_far <- (base_width * 5) - 1

#   i <- 1
#   mean_values <- vector()
#   mean_pos <- vector()
#   while (window_near <= max(pos)) {
#     index <- which(pos > window_near & pos < window_far)
#     if (length(index) >= 10) {
#       value_set <- value[index]
#       pos_set <- pos[index]
#       window_width <- 4
#       if (ave) {
#         for (j in 1:(length(pos_set) - window_width)) {
#           mean_pos[i] <- mean(pos_set[j : (j + window_width)], na.rm = T)
#           mean_values[i] <- mean(value_set[j : (j + window_width)], na.rm = T)
#           i <- i + 1
#         }
#       } else {
#         for (j in 1:(length(pos_set) - window_width)) {
#           mean_pos[i] <- sd(pos_set[j : (j + window_width)], na.rm = T)
#           mean_values[i] <- sd(value_set[j : (j + window_width)], na.rm = T)
#           i <- i + 1
#         }
#       }
#     } else {
#       mean_values[i] <- 0
#       mean_pos[i] <- mean(c(window_near, window_far))
#       i <- i + 1
#     }
#     window_near <- window_near + base_width
#     window_far <- window_far + base_width
#   }

#   for (j in 1:length(mean_values)) {
#     if (is.na(mean_values[j])) {
#       mean_values[j] <- 0
#     }
#   }

#   res <- cbind(mean_pos, mean_values)
#   res <- res[order(res[,1]),]
#   return(res)
# }

EH <- function (genotypes) {
  apply(genotypes, 1, function (row) {
    2 * ((sum(row == 0) / sum(row == 0 | 2)) *
    (sum(row == 2) / sum(row == 0 | 2)))
  })
}