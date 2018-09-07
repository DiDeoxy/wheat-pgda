legends <- function (region, pop.code, colours, title, cex = 0.45, ncol = 1, bg = "white", horiz = F) {
  legend(region, legend = levels(pop.code), pch = 19, col = colours, cex = cex, bg = bg, title = title, ncol = ncol, horiz = horiz)
}

drawRects <- function (pop.code, colourSubset, labelOrder, border) {
  rects <- function(x, y) {
    circos.rect(0.05:(length(labelOrder)-0.95), 
                rep(0,length(labelOrder)), 
                0.95:(length(labelOrder)-0.05), 
                rep(1,length(labelOrder)), 
                border = border, 
                lwd = 0.8,
                col = sapply(as.numeric(pop.code)[labelOrder], function(x) { return(colourSubset[x]) } ))
  }
  
  circos.track(ylim = c(0, 1), 
               panel.fun = rects(x, y), 
               track.height = 0.04, 
               bg.border = NA)
}


