## add in xaxt for last three graphs
source("Analysis\\R\\functions\\funcs_calc_stats.R")
source("Analysis\\R\\functions\\colour_sets.R")

plots2 <- function (positions, values, mainGenes, resiGenes, segments, ab = "both", ylim = c(0, 1), raw = TRUE, smooth = TRUE, plotGenes = TRUE, ave = TRUE, colour = FALSE) {
  par(mfrow = c(7,3), oma = c(4,1,3,3))
  cex <- 0.8
  labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
  count <- 1
  by(cbind(positions, values, mainGenes, resiGenes), segments, function (chr) {
    par(mar = c(0, 2, 0, 0))
    if (count %in% c(19, 20, 21)) { # plot xaxis on last three graphs
      plotter2(chr, count, ylim, raw, smooth, plotGenes, ave, ab, colour = colour, xaxt = "s")
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    } else {
      plotter2(chr, count, ylim, raw, smooth, plotGenes, ave, ab, colour = colour)
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    }
    count <<- count + 1
  })
}

plotter2 <- function (chr, count, ylim = c(0, 1), raw = TRUE, smooth = TRUE, genes = TRUE, ave = TRUE, ab = "both", yaxt = "n", xaxt = "n", selected = FALSE, colour = FALSE) {
  pch <- 20
  colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
  order <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3))
  
  if (raw) {
    if (count %in% seq(3, 21, 3) & yaxt == "n") {
      plot(chr[,1], chr[,2], pch = pch, ylim = ylim, cex = 0.5,
           xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
           col = coloursChroms[order[count]], yaxt = yaxt, xaxt = xaxt)
      axis(4, cex.axis = 0.75, las = 1)
    } else {
      plot(chr[,1], chr[,2], pch = pch, ylim = ylim, cex = 0.5,
           xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
           col = coloursChroms[order[count]], yaxt = yaxt, xaxt = xaxt)
    }
  } else {
    if (count %in% seq(3, 21, 3) & yaxt == "n") {
      plot(chr[,1], chr[,2], pch = pch, ylim = ylim, cex = 0.5,
           xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
           col = "#FFFFFF", yaxt = yaxt, xaxt = xaxt)
      axis(4, cex.axis = 0.75, las = 1)
    } else {
      plot(chr[,1], chr[,2], pch = pch, ylim = ylim, cex = 0.5,
           xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
           col = "#FFFFFF", yaxt = yaxt, xaxt = xaxt)
    }
  }
  
  if (ab == "both") {
    abline(h = signif, col = "red")
    abline(h = signif2, col = "blue")
  } else if (ab == "raw") {
    abline(h = signif, col = "red")
  } else if (ab == "smooth") {
    abline(h = signif2, col = "blue")
  }
 
  if (smooth & colour) {
    smoothed <- sw_calc(chr[,1], chr[,2], ave)
    xx <- c(smoothed[,1], rev(smoothed[,1]))
    yy <- c(rep(0, length(smoothed[,2])), rev(smoothed[,2]))
    polygon(x = xx, y = yy, col = adjustcolor(coloursChroms[order[count]], alpha.f = 0.75), lwd = 2, border = NA)
  } else if (smooth) {
    smoothed <- sw_calc(chr[,1], chr[,2], ave)
    xx <- c(smoothed[,1], rev(smoothed[,1]))
    yy <- c(rep(0, length(smoothed[,2])), rev(smoothed[,2]))
    polygon(x = xx, y = yy, col = adjustcolor("#808080", alpha.f = 0.6), lwd = 2, border = NA)
  }
  
  textCex = 0.8
  if (genes & max(ylim) == 1) {
    # abline(v = chr[,4], col = colourSet[15])
    points(chr[,4], y = rep(0, length(chr[,4])), col = colourSet[15], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,4], y = c(0.1,0.3,0.5,0.7,0.9), labels = chr[,5], pos = 4, srt = 0, col = colourSet[15], cex = textCex)
    } else {
      text(chr[,4], y = c(0.1, 0.1, 0.1, 0.1, 0.1), labels = chr[,5], pos = 4, srt = 90, col = colourSet[15], cex = textCex)
    }
    # abline(v = chr[,7], col = colourSet[19])
    points(chr[,7], y = rep(0, length(chr[,7])), col = colourSet[19], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,7], y = c(0.1,0.3,0.5,0.7,0.9), labels = chr[,8], pos = 4, srt = 0, col = colourSet[19], cex = textCex)
    } else {
      text(chr[,7], y = c(0.1, 0.1, 0.1, 0.1, 0.1), labels = chr[,8], pos = 4, srt = 90, col = colourSet[19], cex = textCex)
    }
  } else {
    # abline(v = chr[,4], col = colourSet[15])
    points(chr[,4], y = rep(0, length(chr[,4])), col = colourSet[15], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,4], y = c(0.05,0.15,0.25,0.35,0.45), labels = chr[,5], pos = 4, srt = 0, col = colourSet[15], cex = textCex)
    } else {
      text(chr[,4], y = c(0.05, 0.05, 0.05, 0.05, 0.05), labels = chr[,5], pos = 4, srt = 90, col = colourSet[15], cex = textCex)
    }
    # abline(v = chr[,7], col = colourSet[19])
    points(chr[,7], y = rep(0, length(chr[,7])), col = colourSet[19], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,7], y = c(0.05,0.15,0.25,0.35,0.45), labels = chr[,8], pos = 4, srt = 0, col = colourSet[19], cex = textCex)
    } else {
      text(chr[,7], y = c(0.05, 0.05, 0.05, 0.05, 0.05), labels = chr[,8], pos = 4, srt = 90, col = colourSet[19], cex = textCex)
    }
  }
}

