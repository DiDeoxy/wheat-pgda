library(extrafont)

read.csv("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\phys_map_dups.csv")
cex <- 0.8
pch <- 20
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
order <- c(rep(1,6), rep(2,6), rep(3,6), rep(4,6), rep(5,6), rep(6,6), rep(7,6))
labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
count <- 1


png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\loci\\duplicates.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
par(mfrow = c(7,6), oma = c(3,3,4,2))
by(duplicates$Position, duplicates$Contig, function (x) {
  if (count %in% seq(1, 5, 2)) { # the first graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot(x, xaxt = "n", yaxt = "n")
    mtext(text = "Long", 3, line = 0.5, cex = cex)
    mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
  } else if (count %in% seq(2, 6, 2)) { # the first graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    if (count == 6) {
      plot(0,xaxt='n',yaxt='n',pch='',ylab='',xlab='', ylim = c(0,1))
      axis(4, at=seq(0, 1, 0.1))
      mtext(text = "Short", 3, line = 0.5, cex = cex)
      par(mar = c(0, 2, 0, 0))
      plot(x, xaxt = "n", yaxt = "n")
      count <<- count + 1
      mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    } else {
      plot(x, xaxt = "n", yaxt = "n")
      axis(4, at=seq(0, 1, 0.1))
      mtext(text = "Short", 3, line = 0.5, cex = cex)
    }
  } else if (count %in% seq(7, 41, 2)) { # the remaining graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    if (count == 17) {
      plot(0,xaxt='n',yaxt='n',pch='',ylab='',xlab='', ylim = c(0,1))
      mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
      par(mar = c(0, 0, 0, 2))
      plot(x, xaxt = "n", yaxt = "n")
      count <<- count + 1
      if (count %in% c(14, 16, 18, 26, 28, 30, 38, 40, 42)) axis(4, at=seq(0, 1, 0.1))
    } else {
      plot(x, xaxt = "n", yaxt = "n")
      mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    }
  } else { # the remaining graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    if (count %in% c(12, 24, 42)) {
      plot(0,xaxt='n',yaxt='n',pch='',ylab='',xlab='', ylim = c(0,1))
      if (count %in% c(14, 16, 18, 26, 28, 30, 38, 40, 42)) axis(4, at=seq(0, 1, 0.1))
      par(mar = c(0, 2, 0, 0))
      plot(x, xaxt = "n", yaxt = "n")
      count <<- count + 1
      mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    } else {
      plot(x, xaxt = "n", yaxt = "n")
      if (count %in% c(14, 16, 18, 26, 28, 30, 38, 40, 42)) axis(4, at=seq(0, 1, 0.1))
    }
  }
  count <<- count + 1
})
dev.off()
