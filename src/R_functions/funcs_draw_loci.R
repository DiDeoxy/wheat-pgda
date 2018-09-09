## add in xaxt for last three graphs
source("Analysis\\R\\R_functions\\funcs_calc_stats.R")
source("Analysis\\R\\R_functions\\colour_sets.R")

plots <- function (positions, values, mainGenes, resiGenes, segments, ylim = c(0, 1),
                   oma = c(4,1,3,3), plotGenes = TRUE, colour = FALSE, ab = FALSE) {
  par(mfrow = c(7,3), oma = oma)
  cex <- 0.8
  labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
  count <- 1
  by(cbind(positions, values, mainGenes, resiGenes), segments, function (chr) {
    par(mar = c(0, 2, 0, 0))
    if (count %in% c(19, 20, 21)) { # plot xaxis on last three graphs
      plotter(chr, count, ylim, plotGenes, colour = colour, xaxt = "s", ab = ab)
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    } else {
      plotter(chr, count, ylim, plotGenes, colour = colour, ab = ab)
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    }
    count <<- count + 1
  })
}

plotter <- function (chr, count, ylim, genes = TRUE, yaxt = "n", xaxt = "n", colour = FALSE, ab = FALSE) {
  pch <- 20
  colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
  order <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3))
  
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
  
  if (ab) {
    abline(v = 7E7)
    abline(v = 3.2E8)
  }
  
  textCex = 0.8
  if (genes) {
    points(chr[,4], y = rep(ylim[1], length(chr[,4])), col = colourSet[15], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,4], y = seq(ylim[1] + 0.01, ylim[2] - 0.01, length.out = 5), labels = chr[,5], pos = 4, srt = 0, col = colourSet[15], cex = textCex)
    } else {
      text(chr[,4], y = rep(ylim[1] + 0.01, 5), labels = chr[,5], pos = 4, srt = 90, col = colourSet[15], cex = textCex)
    }
    points(chr[,7], y = rep(ylim[1], length(chr[,7])), col = colourSet[19], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,7], y = seq(ylim[1] + 0.01, ylim[2] - 0.01, length.out = 5), labels = chr[,8], pos = 4, srt = 0, col = colourSet[19], cex = textCex)
    } else {
      text(chr[,7], y = rep(ylim[1] + 0.01, 5), labels = chr[,8], pos = 4, srt = 90, col = colourSet[19], cex = textCex)
    }
  }
}

plots2 <- function (positions, values, mainGenes, resiGenes, segments, ylim = c(0, 1), plotGenes = TRUE, colour = FALSE) {
  par(mfrow = c(7,3), oma = c(4,1,3,3))
  cex <- 0.8
  labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
  count <- 1
  by(cbind(positions, values, mainGenes, resiGenes), segments, function (chr) {
    
    par(mar = c(0, 2, 0, 0))
    if (count %in% c(19, 20, 21)) { # plot xaxis on last three graphs
      plotter2(chr, count, ylim, plotGenes, colour, xaxt = "s")
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    } else {
      plotter2(chr, count, ylim, plotGenes, colour)
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    }
    count <<- count + 1
  })
}

plotter2 <- function (chr, count, ylim = c(0, 1), genes = TRUE, colour = FALSE, yaxt = "n", xaxt = "n", selected = FALSE) {
  pch <- 20
  colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
  order <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3))
  
  if (count %in% seq(3, 21, 3) & yaxt == "n") {
    plot(chr[,1], chr[,2], pch = pch, ylim = ylim, cex = 0.5,
         xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
         col = "red", yaxt = yaxt, xaxt = xaxt)
    points(chr[,1], y = chr[,3],pch = pch, cex = 0.5, col = "orange")
    axis(4, cex.axis = 0.75, las = 1)
    abline(h = 0, col = "blue", lty = "dotted")
  } else {
    plot(chr[,1], chr[,2], pch = pch, ylim = ylim, cex = 0.5,
         xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
         col = "red", yaxt = yaxt, xaxt = xaxt)
    points(chr[,1], y = chr[,3], pch = pch, cex = 0.5, col = "orange")
    abline(h = 0, col = "blue", lty = "dotted")
  }
  
  textCex = 0.8
  if (genes) {
    points(chr[,5], y = rep(ylim[1], length(chr[,5])), col = colourSet[15], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,5], y = seq(ylim[1] + 0.01, ylim[2] - 0.01, length.out = 5), labels = chr[,6], pos = 4, srt = 0, col = colourSet[15], cex = textCex)
    } else {
      text(chr[,5], y = rep(ylim[1] + 0.01, 5), labels = chr[,6], pos = 4, srt = 90, col = colourSet[15], cex = textCex)
    }
    points(chr[,8], y = rep(ylim[1], length(chr[,5])), col = colourSet[19], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,8], y = seq(ylim[1] + 0.01, ylim[2] - 0.01, length.out = 5), labels = chr[,9], pos = 4, srt = 0, col = colourSet[19], cex = textCex)
    } else {
      text(chr[,8], y = rep(ylim[1] + 0.01, 5), labels = chr[,9], pos = 4, srt = 90, col = colourSet[19], cex = textCex)
    }
  }
}


plots3 <- function (positions, values, mainGenes, resiGenes, segments, ylim = c(0, 1),
                   oma = c(4,1,3,3), plotGenes = TRUE, colour = FALSE, ab = FALSE,
                   haplotypes = NULL) {
  par(mfrow = c(7,3), oma = oma)
  cex <- 0.8
  labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
  count <- 1
  by(cbind(positions, values, mainGenes, resiGenes), segments, function (chr) {
    par(mar = c(0, 2, 0, 0))
    if (count %in% c(19, 20, 21)) { # plot xaxis on last three graphs
      plotter3(chr, count, ylim, plotGenes, colour = colour, xaxt = "s", ab = ab, 
               haplotypes = haplotypes)
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    } else {
      plotter3(chr, count, ylim, plotGenes, colour = colour, ab = ab,
               haplotypes = haplotypes)
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    }
    count <<- count + 1
  })
}

plotter3 <- function (chr, count, ylim, genes = TRUE, yaxt = "n", xaxt = "n",
                      colour = FALSE, ab = FALSE, haplotypes = NULL) {
  pch <- 20
  cex <- 0.5
  colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
  order <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3))
  
  if (count %in% seq(3, 21, 3) & yaxt == "n") {
    plot(chr[,1], chr[,2], pch = pch, ylim = ylim, cex = cex,
         xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
         col = coloursChroms[order[count]], yaxt = yaxt, xaxt = xaxt)
    axis(4, cex.axis = 0.75, las = 1)
  } else {
    colours <- rep(coloursChroms[order[count]], length(chr[,1]))
    plot(chr[,1], chr[,2], pch = pch, ylim = ylim, cex = cex,
         xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
         col = colours, yaxt = yaxt, xaxt = xaxt)
    if (count == 1) {
      points(chr[haplotypes[[1]],1], chr[haplotypes[[1]],2], col = "#000000", pch = 1, cex = 1)
    } else if (count == 4) {
      points(chr[haplotypes[[4]],1], chr[haplotypes[[4]],2], col = "#000000", pch = 1, cex = 1)
    } else if (count == 10) {
      points(chr[haplotypes[[10]],1], chr[haplotypes[[10]],2], col = "#000000", pch = 1, cex = 1)
    } else if (count == 14) {
      points(chr[haplotypes[[14]],1], chr[haplotypes[[14]],2], col = "#000000", pch = 1, cex = 1)
    } else if (count == 16) {
      points(chr[haplotypes[[16]],1], chr[haplotypes[[16]],2], col = "#000000", pch = 1, cex = 1)
    } else if (count == 17) {
      points(chr[haplotypes[[17]],1], chr[haplotypes[[17]],2], col = "#000000", pch = 1, cex = 1)
    } else if (count == 19) {
      points(chr[haplotypes[[19]],1], chr[haplotypes[[19]],2], col = "#000000", pch = 1, cex = 1)
    } 
  }
  
  if (ab) {
    abline(v = 2.5E8)
    abline(v = 3.7E8)
  }
  
  textCex = 0.8
  if (genes) {
    points(chr[,4], y = rep(ylim[1], length(chr[,4])), col = colourSet[15], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,4], y = seq(ylim[1] + 0.01, ylim[2] - 0.01, length.out = 5), labels = chr[,5], pos = 4, srt = 0, col = colourSet[15], cex = textCex)
    } else {
      text(chr[,4], y = rep(ylim[1] + 0.01, 5), labels = chr[,5], pos = 4, srt = 90, col = colourSet[15], cex = textCex)
    }
    points(chr[,7], y = rep(ylim[1], length(chr[,7])), col = colourSet[19], pch = 17, cex = 1)
    if (count %in% c(1,3)) {
      text(chr[,7], y = seq(ylim[1] + 0.01, ylim[2] - 0.01, length.out = 5), labels = chr[,8], pos = 4, srt = 0, col = colourSet[19], cex = textCex)
    } else {
      text(chr[,7], y = rep(ylim[1] + 0.01, 5), labels = chr[,8], pos = 4, srt = 90, col = colourSet[19], cex = textCex)
    }
  }
}