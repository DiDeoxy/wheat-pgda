## add in xaxt for last three graphs
source("Analysis\\R\\R_functions\\funcs_calc_stats.R")
source("Analysis\\R\\R_functions\\colour_sets.R")

plot <- function (positions, values, main_genes, resi_genes, segments,
                   ylim = c(0, 1), oma = c(4, 1, 3, 3), plot_genes = TRUE,
                   colour = FALSE, amova = FALSE, amova_eh = FALSE,
                   all_eh = FALSE) {
  par(mfrow = c(7,3), oma = oma)
  labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"),
                        paste, sep = "")))
  count <- 1
  by(cbind(positions, values, main_genes, resi_genes), segments, 
     function (chr) {
      par(mar = c(0, 2, 0, 0))
      cex <- 0.8
      if (count %in% c(19, 20, 21)) { # plot xaxis on last three graphs
        if (amova) {
          plotter_amova(chr, count, ylim, plot_genes, colour = colour, 
                        xaxt = "s")
        } else if (amova_eh) {
          plotter_amova_eh(chr, count, ylim, plot_genes, colour = colour, 
                           xaxt = "s")
        } else if (all_eh) {
          plotter_all_eh(chr, count, ylim, plot_genes, colour = colour, 
                          xaxt = "s")
        }
        mtext(text = labels[count], 2, line = 0.5, cex = cex)
      } else {
        if (amova) {
          plotter_amova(chr, count, ylim, plot_genes, colour = colour)
        } else if (amova_eh) {
          plotter_amova_eh(chr, count, ylim, plot_genes, colour = colour)
        } else if (all_eh) {
          plotter_all_eh(chr, count, ylim, plot_genes, colour = colour)
        }
        mtext(text = labels[count], 2, line = 0.5, cex = cex)
      }
      count <<- count + 1
    }
  )
}

plotter_amova <- function (chr, count, ylim, genes = TRUE, yaxt = "n", 
                           xaxt = "n", colour = FALSE) {
  pch <- 20
  order <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3),
             rep(5,3), rep(6,3), rep(7,3))

  if (count %in% seq(3, 21, 3) & yaxt == "n") {
    plot(chr[, 1], chr[, 2], pch = pch, ylim = ylim, cex = 0.5,
         xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
         col = colours_chroms[order[count]], yaxt = yaxt, xaxt = xaxt)
    axis(4, cex.axis = 0.75, las = 1)
  } else {
    plot(chr[, 1], chr[, 2], pch = pch, ylim = ylim, cex = 0.5,
         xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
         col = colours_chroms[order[count]], yaxt = yaxt, xaxt = xaxt)
  }
  
  if (genes) {
    genes(chr, ylim)
  }
}

genes <- function(chr, ylim) {
  text_cex <- 0.8
  points(chr[,4], y = rep(ylim[1], length(chr[,4])), col = colour_set[15],
         pch = 17, cex = 1)
  if (count %in% c(1,3)) {
    text(chr[,4], y = seq(ylim[1] + 0.01, ylim[2] - 0.01, length.out = 5),
         labels = chr[,5], pos = 4, srt = 0, col = colour_set[15],
         cex = text_cex)
  } else {
    text(chr[,4], y = rep(ylim[1] + 0.01, 5), labels = chr[,5], pos = 4,
         srt = 90, col = colour_set[15], cex = text_cex)
  }
  points(chr[,7], y = rep(ylim[1], length(chr[,7])), col = colour_set[19],
         pch = 17, cex = 1)
  if (count %in% c(1,3)) {
    text(chr[,7], y = seq(ylim[1] + 0.01, ylim[2] - 0.01, length.out = 5),
         labels = chr[,8], pos = 4, srt = 0, col = colour_set[19],
         cex = text_cex)
  } else {
    text(chr[,7], y = rep(ylim[1] + 0.01, 5), labels = chr[,8], pos = 4,
         srt = 90, col = colour_set[19], cex = text_cex)
  }
}

plotter_amova_eh <- function (chr, count, ylim = c(0, 1), genes = TRUE, 
                       colour = FALSE, yaxt = "n", xaxt = "n", 
                       selected = FALSE) {
  pch <- 20
  order <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3),
             rep(5,3), rep(6,3), rep(7,3))
  
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
  
  if (genes) {
    genes(chr, ylim)
  }
}

plotter_all_eh <- function (chr, count, ylim, genes = TRUE, yaxt = "n", 
                            xaxt = "n", colour = FALSE, ab = FALSE,
                            haplotypes = NULL) {
  pch <- 20
  cex <- 0.5
  if (count %in% seq(3, 21, 3) & yaxt == "n") {
    plot(chr[, 1], chr[, 2], pch = pch, ylim = ylim, cex = cex,
         xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
         col = colours_chroms[order[count]], yaxt = yaxt, xaxt = xaxt)
    axis(4, cex.axis = 0.75, las = 1)
  } else {
    colours <- rep(colours_chroms[order[count]], length(chr[, 1]))
    plot(chr[, 1], chr[, 2], pch = pch, ylim = ylim, cex = cex,
         xlim = c(0, max(snp.pos)), cex.axis = 0.8, xlab = "", ylab = "",
         col = colours, yaxt = yaxt, xaxt = xaxt)
    if (count == 1) {
      points(chr[haplotypes[[1]],1], chr[haplotypes[[1]],2], col = "#000000",
             pch = 1, cex = 1)
    } else if (count == 4) {
      points(chr[haplotypes[[4]],1], chr[haplotypes[[4]],2], col = "#000000",
             pch = 1, cex = 1)
    } else if (count == 10) {
      points(chr[haplotypes[[10]],1], chr[haplotypes[[10]],2], col = "#000000",
             pch = 1, cex = 1)
    } else if (count == 14) {
      points(chr[haplotypes[[14]],1], chr[haplotypes[[14]],2], col = "#000000",
             pch = 1, cex = 1)
    } else if (count == 16) {
      points(chr[haplotypes[[16]],1], chr[haplotypes[[16]],2], col = "#000000",
             pch = 1, cex = 1)
    } else if (count == 17) {
      points(chr[haplotypes[[17]],1], chr[haplotypes[[17]],2], col = "#000000",
             pch = 1, cex = 1)
    } else if (count == 19) {
      points(chr[haplotypes[[19]],1], chr[haplotypes[[19]],2], col = "#000000",
             pch = 1, cex = 1)
    } 
  }

  if (genes) {
    genes(chr, ylim)
  }
}