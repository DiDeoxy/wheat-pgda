#' Output figures based on map stats
#'
#' Outputs a figure of gap distances and neighbouring ld by genome and overall
#'
#' @param lng data supplied by the calc_lng function
#' @param genome_ld data supplied by the calc_cld_stats function
#' @param gds the gds object
#' @param plot_title the title of the plot
#' @param y_lim the y-axis limit of the plot
#' @param out_name the name for the output file
#'
#' @importFrom dplyr tibble
#' @importFrom GGally ggmatrix
#' @importFrom ggplot2 aes element_text geom_freqpoly ggplot scale_colour_manual theme xlim ylim
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_c
#'
#' @return None

plot_gaps_nbs_ld <- function(lng, genome_ld, gds, plot_title, y_lim, out_name) {
  # histograms and boxplots depicting the distribution of gaps on each genome
  gaps_log10 <- tibble(
    Genome = factor(
      c(rep("A", length(lng$A$gaps)),
        rep("B", length(lng$B$gaps)),
        rep("D", length(lng$D$gaps)),
        rep("All", length(c(lng$A$gaps, lng$B$gaps, lng$D$gaps)))
      ),
      levels = c("A", "B", "D", "All")
    ),
    gaps = c(
      lng$A$gaps,
      lng$B$gaps,
      lng$D$gaps,
      lng$A$gaps, lng$B$gaps, lng$D$gaps
    ) %>% log10()
  )

  # histograms and boxplots depicting the distribution of gaps on each genome
  nbs_ld_genome <- tibble(
    Genome = factor(
      c(rep("A", length(genome_ld$A$nbs)),
        rep("B", length(genome_ld$B$nbs)),
        rep("D", length(genome_ld$D$nbs)),
        rep("All", length(c(genome_ld$A$nbs, genome_ld$B$nbs, genome_ld$D$nbs)))
      ),
      levels = c("A", "B", "D", "All")
    ),
    ld = c(
      genome_ld$A$nbs,
      genome_ld$B$nbs,
      genome_ld$D$nbs,
      genome_ld$A$nbs, genome_ld$B$nbs, genome_ld$D$nbs
    )
  )

  plots <- list()
  plots[[1]] <- gaps_log10 %>%
    ggplot() +
    geom_freqpoly(aes(gaps, colour = Genome), size = 0.3) +
    scale_colour_manual(values = brewer.pal(4, "Dark2")) +
    xlim(min(gaps_log10$gaps), max(gaps_log10$gaps)) +
    ylim(0, y_lim)
  plots[[2]] <- nbs_ld_genome %>%
    ggplot() +
    geom_freqpoly(aes(ld, colour = Genome), size = 0.3) +
    scale_colour_manual(values = brewer.pal(4, "Dark2")) +
    xlim(0, 1.0001) +
    ylim(0, y_lim)

  # turn plot list into ggmatrix
  plots_matrix <- ggmatrix(
    plots, nrow = 1, ncol = 2, yAxisLabels = "Num Markers",
    xAxisLabels = c(
      "Log 10 of Mb Gap Distances",
      "Abs. LD Between Neighbouring Markers"
    ),
    title = plot_title,
    legend = c(1, 2)
  )

  # plot the matrix
  png(
    file.path(
      map_stats_and_plots, str_c(basename(gds), ".png")
    ),
    family = "Times New Roman", width = 100, height = 62, pointsize = 10,
    units = "mm", res = 300)
  print(plots_matrix + 
    theme(legend.position = "bottom",
      text = element_text(
        size = 8, lineheight = 0.1)))
  dev.off()
}