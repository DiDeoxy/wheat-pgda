library(adegenet)
library(poppr)
library(tidyverse)

i <- 20
j <- 0.1

# phis <- list()
# for (i in c(20, 40, 80, 160, 320, 640)) {
#   half <- i / 2
#   genotypes <- tibble(g = c(rep("A", half), rep("B", half)))
#   for (j in seq(0, 0.5, 0.1)) {
#     strata <- tibble(s = 
#     c(rep("A", half * j), rep("B", half - (half * j)),
#       rep("B", half * j), rep("A", half - (half * j))
#       # rep("A", half)
#       )
#     )

#     test_genind <- df2genind(
#       genotypes, ind.names = as.character(1:i), loc.names = "x",
#       ploidy = 1, type = "codom", ncode = 1, strata = strata
#     )

#     if (length(phis[[as.character(i)]])) {
#       phis[[as.character(i)]] <- c(
#         phis[[as.character(i)]],
#         poppr.amova(test_genind, hier = ~s,
#           missing = "genotype", within = F, clonecorrect = F
#         )$statphi
#       )
#     } else {
#       phis[[as.character(i)]] <- poppr.amova(
#         test_genind, hier = ~s, missing = "genotype", within = F,
#         clonecorrect = F
#       )$statphi
#     }
#   }
#   phis[[as.character(i)]] <- rev(phis[[as.character(i)]])
# }
# row_top <- c(rep(0, 182), rep(1, 182))
# row_bot <- c(rep(1, 182), rep(0, 182))
# sample_dist <- matrix(c(rep(row_top, 182), rep(row_bot, 182)), nrow = 364) %>%
#   as.dist()

# long_phis <- unlist(phis) %>% matrix(ncol = 6) %>% as.tibble() %>%
#   select("20" = V1, "40" = V2, "80" = V3, "160" = V4, "320" = V5,
#     "640" = V6) %>%
#   gather("Sample Size", phi, 
#     c("20", "40", "80", "160", "320", "640"))

# long_phis <- long_phis %>% add_column(Percent = rep(seq(0, 50, 10), 6))
# long_phis["Sample Size"] <- factor(long_phis[["Sample Size"]], 
#   levels = c("20", "40", "80", "160", "320", "640"))

# long_phis %>%
#   ggplot() +
#     geom_line(aes(Percent, phi, colour = !!as.symbol("Sample Size")))

factors <- function(x) {
    x <- as.integer(x)
    div <- seq_len(abs(x))
    factors <- div[x %% div == 0L]
    factors[2:(length(factors) - 1)]
}
factors(26)

num_samples <- 364
half <- num_samples / 2
factors_half <- factors(half)
mid_factor <- factors_half[(length(factors_half) / 2) + 2]
strata <- tibble(s = c(rep("top", half), rep("bot", half)))
phi_data <- tibble(num_A = integer(), num_A_top = integer(), phi = double())
for (num_A in seq(mid_factor, half, mid_factor)) {
  nums_A_top <- 0
  # if ((num_A / 2) %% 2) {
  #   nums_A_top <- seq(1, num_A / 2, 2)
  # } else {
    nums_A_top <- seq(0, num_A / 2, 2)
  # }
  for (num_A_top in nums_A_top) {
    genotypes <- tibble(g =
    c(rep("A", num_A_top), rep("B", half - num_A_top),
      rep("A", num_A - num_A_top), rep("B", half - (num_A - num_A_top))
      )
    )
    # print(nrow(genotypes))

    test_genind <- df2genind(
      genotypes, ind.names = as.character(1:num_samples), loc.names = "x",
      ploidy = 1, type = "codom", ncode = 1, strata = strata
    )
    # print(test_genind)

    phi_data <- phi_data %>% add_row(num_A = num_A, num_A_top = num_A_top,
      phi = poppr.amova(
        test_genind, hier = ~s, missing = "genotype", clonecorrect = F,
        within = F
      )$statphi[[1]]
    )
  }
}
phi_data

phi_data$num_A <- factor(phi_data$num_A)

phi_data %>%
  ggplot() +
    geom_line(aes(num_A_top, phi, colour = num_A))