#' Craws rectangles for circliz graph
#'
#' Function for drawing rectangles in each track of a cirlcize graph
#'
#' @importFrom circlize circos.rect circos.track
#'
#' @param pop_code a factor of the classifications of the individuals
#' @param colour_subset the colours used for the factor levels
#' @param the order of the factors in the row
#' @param border the colour for the border of the rectangles
#'
#' @return None
#'
#' @export

draw_rects <- function(pop_code, colour_subset, label_order, border) {
  leng <- length(label_order)
  rects <- function(x, y) {
    circos.rect(
      0.05:(leng - 0.95), rep(0, leng),
      0.95:(leng - 0.05), rep(1, leng),
      border = border,
      lwd = 0.9,
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