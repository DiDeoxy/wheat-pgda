# genind object for amova
strata_genind <- read_rds(file.path(geninds, "strata.rds"))

# create an IBS distance object of the sample genotypes
wheat_gds <- snpgdsOpen(ld_gds)
sample_dist <- as.dist(1 - snpgdsIBS(wheat_gds, autosome.only = F)$ibs)
snpgdsClose(wheat_gds)

strata <- c(
  "bp", "era", "pheno", "pheno/bp", "pheno/era", "pheno/bp/era", "pheno/era/bp", "bp_era",
  "pheno/bp_era", "clusters", "clusters/bp_era"
)

amova_table <- tibble(
  Comparison = character(), `% Variation` = double(), `p-value` = double()
)
for (stratum in strata) {
  print(stratum)
  amova_result <- poppr.amova(
    strata_genind, hier = str_c("~", stratum) %>% as.formula(), cutoff = 0.1,
    missing = "genotype", within = FALSE, clonecorrect = FALSE, 
    dist = sample_dist, quiet = T
  )
  amova_randtest <- NULL
  if (nrow(amova_result$statphi) == 1) {
    amova_randtest <- randtest(amova_result, nrepet = 999, alter = "greater")
    amova_table <- amova_table %>%
      add_row(
        Comparison = stratum,
        `% Variation` = round(amova_result$componentsofcovariance$`%`[1], 2),
        `p-value` = amova_randtest$pvalue
      )
  } else if (nrow(amova_result$statphi) == 3) {
    amova_randtest <- randtest(amova_result, nrepet = 999, output = "full") %>%
    with(
      as.krandtest(
        sim, obs, alter = c("greater", "greater"), call = call,
        names = names
      )
    )
    amova_table <- amova_table %>%
      add_row(
        Comparison = str_c(stratum, " (", strsplit(stratum, "/")[[1]][1], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[1], 2),
        `p-value` = amova_randtest$pvalue[3]
      )
    amova_table <- amova_table %>%
      add_row(
        Comparison = str_c(stratum, " (", strsplit(stratum, "/")[[1]][2], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[2], 2),
        `p-value` = amova_randtest$pvalue[2]
      )
  } else {
    amova_randtest <- randtest(amova_result, nrepet = 999, output = "full") %>%
    with(
      as.krandtest(
        sim, obs, alter = c("greater", "greater", "greater"), call = call,
        names = names
      )
    )
    amova_table <- amova_table %>%
      add_row(
        Comparison = str_c(stratum, " (", strsplit(stratum, "/")[[1]][1], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[1], 2),
        `p-value` = amova_randtest$pvalue[3]
      )
    amova_table <- amova_table %>%
      add_row(
        Comparison = str_c(stratum, " (", strsplit(stratum, "/")[[1]][2], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[2], 2),
        `p-value` = amova_randtest$pvalue[2]
      )
    amova_table <- amova_table %>%
      add_row(
        Comparison = str_c(stratum, " (", strsplit(stratum, "/")[[1]][3], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[3], 2),
        `p-value` = amova_randtest$pvalue[2]
      )
  }
}
amova_table

write_csv(
  amova_table, file.path("results", "amova_table.csv"),col_names = TRUE
)

# # correlation stuff
# mis_pheno <- which(strata_genind$strata$pheno == "UNKNOWN")
# mis_bp <- which(strata_genind$strata$bp == "N/A")
# mis_era <- which(strata_genind$strata$era == "UNKNOWN")
# mis_pheno_bp <- unique(c(mis_pheno, mis_bp))
# mis_pheno_era <- unique(c(mis_pheno, mis_era))
# mis_bp_era <- unique(c(mis_era, mis_bp))
# mis_pheno_bp_era <- unique(c(mis_bp_era, mis_pheno))
# table(strata_genind$strata$bp[-mis_pheno_bp], strata_genind$strata$pheno[-mis_pheno_bp])
# cor.test(strata_genind$strata$bp[-mis_pheno_bp] %>% as.integer(), strata_genind$strata$pheno[-mis_pheno_bp] %>% as.integer())
# table(strata_genind$strata$era[-mis_pheno_era], strata_genind$strata$pheno[-mis_pheno_era])
# cor.test(strata_genind$strata$era[-mis_pheno_era] %>% as.integer(), strata_genind$strata$pheno[-mis_pheno_era] %>% as.integer())
# table(strata_genind$strata$bp[-mis_bp_era], strata_genind$strata$era[-mis_bp_era])
# cor.test(strata_genind$strata$bp[-mis_bp_era] %>% as.integer(), strata_genind$strata$era[-mis_bp_era] %>% as.integer())
# table(strata_genind$strata$bp_era[-mis_pheno_bp_era], strata_genind$strata$pheno[-mis_pheno_bp_era])
# cor.test(strata_genind$strata$bp[-mis_pheno_bp_era] %>% as.integer(), strata_genind$strata$era[-mis_pheno_bp_era] %>% as.integer())