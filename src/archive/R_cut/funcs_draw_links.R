createLinkPlot <- function (folder, sq, bandWidth, group1, group2, plotDesig = FALSE) {
  ## circos diagram
  for (base in sq) {
    if (group2 == "all") {
      png(paste("Results\\links\\", folder, "\\",
                "between_", group1, "_", sample.id[labelOrder][group1], "_and_", group2,
                "_at_freq_", as.character(i), "_to_", as.character(i+sqPlus),
                ".png", sep = ""),
          family="Times New Roman", width = 11, height = 11, pointsize = 20, units = "in", res = 500)
    } else {
      png(paste("Results\\links\\", folder, "\\",
                "between_", group1, "_and_", group2,
                "_at_freq_", as.character(i), "_to_", as.character(i+sqPlus),
                ".png", sep = ""),
          family="Times New Roman", width = 11, height = 11, pointsize = 20, units = "in", res = 500)
    }

    circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0, track.margin = c(0.005, 0.005))
    circos.initialize(as.character(1:358), xlim = c(0, 1))
    circos.track(ylim = c(0, 1), track.height = 0.15, bg.border = NA)
    circos.track(ylim = c(0, 1), track.height = 0.04, bg.border = NA)

    for (j in 1:358) {
      circos.text(1, 0, labels = sample.id[labelOrder][j], sector.index = j, track.index = 1,
                  facing = "clockwise", niceFacing = T, cex = 0.25, adj = c(0, -0.2), font = 2)
    }

    if (plotDesig == TRUE) {
      circos.track(ylim = c(0, 1), track.height = 0.04, bg.border = NA)
      for (j in 1:358) {
        circos.rect(0.05, 0, 0.95, 1, sector.index = j, track.index = 2, border = "#808080",
                    lwd = 0.8, col = coloursDesig[as.numeric(desig)][labelOrder][j])
      }
      for (j in 1:358) {
        circos.rect(0.05, 0, 0.95, 1, sector.index = j, track.index = 3, border = "#808080",
                    lwd = 0.8, col = coloursIcl[as.numeric(bests)][labelOrder][j])

      }
    } else {
      for (j in 1:358) {
        circos.rect(0.05, 0, 0.95, 1, sector.index = j, track.index = 2, border = "#808080",
                    lwd = 0.8, col = coloursIcl[as.numeric(bests)][labelOrder][j])

      }
    }
    print(c(i, i+sqPlus))
    topIndices <- which(totalDistances >= quantile(totalDist, probs = base, na.rm = T)
                        & totalDistances <= quantile(totalDist, probs = base+bandWidth, na.rm = T), arr.ind = T)

    links(topIndices, group1, group2)
    circos.clear()
    dev.off()
  }
}

links <- function (topIndices, group1, group2) {
  apply(topIndices, 1, function (row) {
    if (group2 == "all") {
      if (row[1] == group1 |
          row[2] == group1) {
        circos.link(row[1], 0.5, row[2], 0.5)
      }
    } else {
      if (bests[labelOrder][row[1]] == group1 & bests[labelOrder][row[2]] == group2 |
          bests[labelOrder][row[1]] == group2 & bests[labelOrder][row[2]] == group1) {
        circos.link(row[1], 0.5, row[2], 0.5)
      }
    }
  })
}

specificLinks <- function (base, bandWidth, indiv) {
  topIndices <- which(totalDistances >= quantile(totalDist, probs = base, na.rm = T) &
                      totalDistances <= quantile(totalDist, probs = base+bandWidth, na.rm = T), arr.ind = T)
  links <- vector()
  for (row in 1:nrow(topIndices))  {
    if (topIndices[row,1] == indiv) {
      links <- c(links, topIndices[row,2])
    } else if (topIndices[row,2] == indiv) {
      links <- c(links, topIndices[row,1])
    }
  }
  return(links)
}

averageLinks <- function (base, bandWidth, indiv) {
  topIndices <- which(totalDistances >= quantile(totalDist, probs = base, na.rm = T) &
                        totalDistances <= quantile(totalDist, probs = base+bandWidth, na.rm = T), arr.ind = T)
  links <- vector()
  for (row in 1:nrow(topIndices))  {
    if (topIndices[row,1] == indiv) {
      links <- c(links, percentIBD[topIndices[row,][1], topIndices[row,][2]])
    } else if (topIndices[row,2] == indiv) {
      links <- c(links, percentIBD[topIndices[row,][1], topIndices[row,][2]])
    }
  }
  if (length(links) > 10) {
    return(sum(links)/length(links))
  } else {
    return(NA)
  }
  
}


potentialConnections <- function(dim1, dim2) {
  connectionsKept <- apply(expand.grid(1:dim1, 1:dim2), 1, function (row) {
    if (row[1] >= row[2]) {
      return(row)
    }
  })
  unlistedConnectionsKept <- data.frame()
  for (row in connectionsKept) {
    if (!is.null(row)) {
      unlistedConnectionsKept <- rbind(unlistedConnectionsKept, row)
    }
  }
  unlistedConnectionsKept
}

# potentialConnections(4, 3)

linkage <- function (floor, ceiling, end1, end2) {
  topIndices <- which(totalDistances >= quantile(totalDist, probs = floor, na.rm = T) &
                      totalDistances <= quantile(totalDist, probs = ceiling, na.rm = T), 
                      arr.ind = T)
  linkPresence <- apply(topIndices, 1, function (row) {
    if (bests[labelOrder][row[1]] == end1 & bests[labelOrder][row[2]] == end2 |
        bests[labelOrder][row[1]] == end2 & bests[labelOrder][row[2]] == end1) {
      return(1)
    } else {
      return(0)
    }
  })
  # print(sum(linkPresence))
  # potCons <- nrow(potentialConnections(length(which(bests == end1)), length(which(bests == end2))))
  numPerms <- length(which(bests == end1)) * length(which(bests == end2))
  # print(potCons)
  if (end1 == end2) {
    return(paste("Percent of possbile connections between ", end1, " and ", end2, " from ", floor, " to ", ceiling, " ", round((sum(linkPresence)/numPerms)*100, 2)))
  } else {
    return(paste("Percent of possbile connections between ", end1, " and ", end2, " from ", floor, " to ", ceiling, " ", round((sum(linkPresence)/numPerms)*50, 2)))
  }
}

# linkage(0, 1, 1, 2)

# hist(as.vector(totalDistances))
# which(totalDistances >= quantile(totalDist, probs = i, na.rm = T) & totalDistances <= quantile(totalDist, probs = i+0.005, na.rm = T), arr.ind = T)
# circos.track(ylim = c(0, 1), track.height = 0.04, bg.border = NA)
# circos.track(ylim = c(0, 1), track.height = 0.04, bg.border = NA)
# circos.track(ylim = c(0, 1), track.height = 0.04, bg.border = NA)
# circos.track(ylim = c(0, 1), track.height = 0.04, bg.border = NA)

# for (i in 1:358) {
#   circos.rect(0.05, 0, 0.95, 1, sector.index = i, track.index = 2, border = "#808080", 
#               lwd = 0.8, col = coloursEra[as.numeric(era)[labelOrder][i]])
# }
# for (i in 1:358) {
#   circos.rect(0.05, 0, 0.95, 1, sector.index = i, track.index = 3, border = "#808080", 
#               lwd = 0.8, col = coloursBp[as.numeric(bp)[labelOrder][i]])
# }
# for (i in 1:358) {
#   circos.rect(0.05, 0, 0.95, 1, sector.index = i, track.index = 4, border = "#808080", 
#               lwd = 0.8, col = coloursMc[as.numeric(mc)[labelOrder][i]])
# }
# for (i in 1:358) {
#   circos.rect(0.05, 0, 0.95, 1, sector.index = i, track.index = 5, border = "#808080", 
#               lwd = 0.8, col = coloursDesig[as.numeric(desig)[labelOrder][i]])
# }