library(ggplot2)
library(plotly)
library(Rtsne)
import::from(ape, "pcoa")
import::from(dplyr, "bind_cols")
import::from(htmlwidgets, "saveWidget")
import::from(purrr, "map", "map_dfc", "pmap", "reduce")
import::from(magrittr, "%>%")
import::from(readr, "read_csv")
import::from(stringr, "str_c")
import::from(tibble, "tibble")

# simulate major allele frequencies
# set.seed(100)
freqs <- runif(6, 0, 1)
mjafs <- pmax(freqs, 1 - freqs)

# calculate genotype frequencies
gfs <- tibble(
  p = mjafs ^ 2,
  pq = 2 * (mjafs * (1 - mjafs)),
  q = (1 - mjafs) ^ 2
)

pop <- pmap(gfs, list) %>% expand.grid() %>% data.matrix()

geno_freq <- pmap(pop %>% as.data.frame(), prod) %>% as.numeric()

gvs <- map_dfc(seq_along(mjafs), function (i) {
  sample.int(3, 3)
}) %>% t()

popvs <- pmap(gvs %>% as.data.frame(), list) %>% expand.grid() %>% data.matrix()

popvs_mean <- rowMeans(popvs)

pd <- dist(pop)

gfsd <- geno_freq %>% sort()


################################################################################

histo <- plot_ly(
  x = log10(geno_freq),
  type = "histogram",
  cumulative = list(enabled = TRUE, direction = "decreasing")
)

saveWidget(
  as_widget(histo),
  "/home/maxh/projects/wheat-pgda/results/cdf.html"
)

histo <- plot_ly(
  x = log10(geno_freq),
  type = "histogram"
)

saveWidget(
  as_widget(histo),
  "/home/maxh/projects/wheat-pgda/results/histo.html"
)

################################################################################

pca <- prcomp(pd)

tsne <- Rtsne(pd, dims = 3, perplexity = 20)

pco <- pcoa(pd)

values <- pca$x %>% as.data.frame()
values2 <- tsne$Y %>% as.data.frame()
values3 <- pco$vectors %>% as.data.frame()

scatter <- plot_ly() %>% 
  add_markers(x = values$PC1, y = values$PC2, z = values$PC3, color = geno_freq)

saveWidget(
  as_widget(scatter),
  "/home/maxh/projects/wheat-pgda/results/scatter_pca.html"
)

scatter <- plot_ly() %>% 
  add_markers(x = values2[,1], y = values2[,2], z = values2[,3], color = geno_freq)

saveWidget(
  as_widget(scatter),
  "/home/maxh/projects/wheat-pgda/results/scatter_tsne.html"
)

scatter <- plot_ly() %>% 
  add_markers(x = values3[,1], y = values3[,2], z = values3[,3], color = geno_freq)

saveWidget(
  as_widget(scatter),
  "/home/maxh/projects/wheat-pgda/results/scatter_pcoa.html"
)

scatter <- plot_ly() %>% 
  add_trace(type = 'histogram2dcontour', x = popvs_mean, y = log10(geno_freq)) %>%
  add_trace(type = "scatter", x = popvs_mean[you], y = log10(geno_freq)[you])

saveWidget(
  as_widget(scatter),
  "/home/maxh/projects/wheat-pgda/results/density2d_values_freq.html"
)