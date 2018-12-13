draw_rects <- function(pop_code, colour_subset, label_order, border) {
  leng <- length(label_order)
  rects <- function(x, y) {
    circos.rect(
      0.05:(leng - 0.95), rep(0, leng),
      0.95:(leng - 0.05), rep(1, leng),
      border = border,
      lwd = 1.2,
      col = sapply(
        as.numeric(pop_code)[label_order],
        function(x) {
          return(colour_subset[x])
        }
      )
    )
  }

  circos.track(
    ylim = c(0, 1),
    panel.fun = rects(x, y),
    track.height = 0.04,
    bg.border = NA
  )
}