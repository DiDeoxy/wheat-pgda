source(file.path("src", "R", "file_paths.R"))
library(poppr)
library(tidyverse)
library(parallel)
# install.packages("matrixStats")
library(matrixStats)
import::from(pgda, "snpgds_parse")

wheat_data <- snpgds_parse(phys_gds)



clusters <- factor(read_rds(hdbscan)$cluster)
levels(clusters) <- c(
  "Noise", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"
)

categorizations <- list(
  clusters, wheat_data$sample$annot$bp, wheat_data$sample$annot$era,
  wheat_data$sample$annot$mc, wheat_data$sample$annot$pheno,
  wheat_data$sample$annot$habit, wheat_data$sample$annot$colour,
  wheat_data$sample$annot$texture
)
categorization_names <- c(
  "HDBSCAN Clusters", "Breeding Program", "Era", "Market Class", "Phenotype",
  "Growth Habit", "Colour", "Texture"
)
categorization_colours <- list(
  colours_hdbscan_pic_legend, colours_bp_pic, colours_era, colours_mc, colours_pheno_pic,
  colour_set[c(1, 22, 4)], colour_set[c(1, 22, 4)], colour_set[c(1, 4, 22)]
)

lapply(1, function (i) {
  categorization <- categorizations[[i]]
  print(categorization_names[i])
  categories <- categorization %>% as.factor() %>% levels()

  rarefied <- lapply(seq_along(categories), function (j) {
    category <- categories[j]
    indivs <- which(categorization == category)
    subset_sizes <- lseq(2, max(2, length(indivs)), 10) %>% floor() %>% unique()
    if (length(subset_sizes) > 1) {
      lapply(subset_sizes, function (subset_size) {
        print(str_c(category, ": ", subset_size))
        tibble(
          category = category, subset_size = subset_size,
          ar = allele.richness(wheat_data$geno[, indivs], subset_size)
        )
      })
    }
  })
})


allele.richness <- function (pop, size) {
  p_q <- cbind(rowSums(pop == 0), rowSums(pop == 2))
  n <- rowSums(p_q)
  mclapply(1:length(n), function (marker) {
    (1 - lapply(0:(size - 1), function(g) {
      (n[marker] - p_q[marker, ] - g) / (n[marker] - g)
    }) %>% do.call(rbind, .) %>% colProds()) %>% sum()
  }, mc.cores = detectCores()) %>% unlist() %>% mean()
}

allele.richness <- function (pop) {
  p_q <- cbind(rowSums(pop == 0), rowSums(pop == 2))
  n <- rowSums(p_q)
  lapply(1:length(n), function (marker) {
    dif <- n[marker] - p_q[marker, ]
    inter <- lapply(0:(ncol(pop) - 1), function (g) {
      if (2 == sum(n[marker] > p_q[marker, ] + g)) {
        (dif - g) / (n[marker] - g)
      }
    }) %>% compact() %>% do.call(rbind, .)
    print(inter)
    # lapply(1:nrow(temp), function (row) {
    #   print(temp[1:row, ])
    #   # (1 - colProds(temp[1:row, ])) %>% sum()
    # })
  })
}

allele.richness(wheat_data$geno[, which(wheat_data$sample$annot$era == "2001-2016")])

, mc.cores = detectCores()
