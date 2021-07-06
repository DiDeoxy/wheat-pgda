source("wheat-pgda/src/R/15_chrom_plots/1_base_data.R")

text_size <- 9

mono_haplos <- cbind(
  id = "mono_haplo", type = "mono_haplo", base = 0.5,
  read_csv(file.path(intermediate, "haplo_windows.csv"))
) %>%
  split(.$chrom)

plots <- lapply(names(mono_haplos), function (chrom) {
  chrom_data <- snp_data[[chrom]] %>%
    rowwise() %>%
    mutate(
      Type = c("All Other", "Long Haplotype")[
        which(
          c(
            pos_mb < mono_haplos[[chrom]]$start |
            pos_mb > mono_haplos[[chrom]]$end,
            pos_mb >= mono_haplos[[chrom]]$start &
            pos_mb <= mono_haplos[[chrom]]$end
          )
        )
      ]
    )

  chrom_data %>% ggplot(aes(Type, eh)) +
    ylim(0, 0.5) +
    geom_violin(aes(fill = Type)) +
    geom_boxplot(width = 0.2) +
    labs(y = "Diversity", title = chrom) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = text_size * 3),
      axis.text.y = element_text(size = text_size * 3),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = text_size * 3),
      plot.title = element_text(size = text_size * 3)
    )
})

png(
  file.path("results", "mono_haplo_vs_distal_eh_dists.png"),
  family = "Times New Roman", width = 1240, height = 3508
)
grid.arrange(
  grobs = plots, nrow = 3, ncol = 2
)
dev.off()