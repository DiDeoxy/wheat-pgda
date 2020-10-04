source("wheat-pgda/src/R/file_paths.R")
import::from(ade4, "as.krandtest", "quasieuclid", "randtest")
import::from(dplyr, "add_row")
import::from(magrittr, "%>%")
import::from(poppr, "poppr.amova")
import::from(readr, "read_rds", "read_csv", "write_csv")
import::from(SNPRelate, "snpgdsClose", "snpgdsIBS", "snpgdsOpen")
import::from(stringr, "str_c")
import::from(tibble, "add_column", "tibble")
import::from(vcd, "assocstats")

# genind object for amova
strata_genind <- read_rds(file.path(geninds, "strata.rds"))

# create an IBS distance object of the sample genotypes
wheat_gds <- snpgdsOpen(ld_phys_gds)
sample_dist <- as.dist(1 - snpgdsIBS(wheat_gds, autosome.only = F)$ibs)
snpgdsClose(wheat_gds)

hierarchies <- c(
  "bp", "era", "mtg", "mtg/bp", "mtg/era", "mtg/bp/era", "mtg/era/bp",
  "bp_era", "mtg/bp_era", "clusters", "clusters/bp_era"
)

amova_table <- tibble(
  Hierarchy = character(), `% Variation` = double(), `p-value` = double()
)
for (hierarchy in hierarchies) {
  print(hierarchy)
  amova_result <- poppr.amova(
    strata_genind, hier = str_c("~", hierarchy) %>% as.formula(), cutoff = 0.1,
    missing = "genotype", within = FALSE, clonecorrect = FALSE, 
    dist = sample_dist, quiet = T
  )
  amova_randtest <- NULL
  if (nrow(amova_result$statphi) == 1) {
    amova_randtest <- randtest(amova_result, nrepet = 999, alter = "greater")
    amova_table <- amova_table %>%
      add_row(
        Hierarchy = hierarchy,
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
        Hierarchy = str_c(hierarchy, " (", strsplit(hierarchy, "/")[[1]][1], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[1], 2),
        `p-value` = amova_randtest$pvalue[3]
      )
    amova_table <- amova_table %>%
      add_row(
        Hierarchy = str_c(hierarchy, " (", strsplit(hierarchy, "/")[[1]][2], ")"),
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
        Hierarchy = str_c(hierarchy, " (", strsplit(hierarchy, "/")[[1]][1], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[1], 2),
        `p-value` = amova_randtest$pvalue[3]
      )
    amova_table <- amova_table %>%
      add_row(
        Hierarchy = str_c(hierarchy, " (", strsplit(hierarchy, "/")[[1]][2], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[2], 2),
        `p-value` = amova_randtest$pvalue[2]
      )
    amova_table <- amova_table %>%
      add_row(
        Hierarchy = str_c(hierarchy, " (", strsplit(hierarchy, "/")[[1]][3], ")"),
        `% Variation` = round(amova_result$componentsofcovariance$`%`[3], 2),
        `p-value` = amova_randtest$pvalue[2]
      )
  }
}
amova_table

amova_table$Hierarchy <- c(
  "BP", "Era", "MTG", "MTG(BP): MTG", "MTG(BP): BP", "MTG(Era): MTG",
  "MTG(Era): Era", "MTG(BP(Era)): MTG", "MTG(BP(Era)): BP", "MTG(BP(Era)): Era",
  "MTG(Era(BP)): MTG", "MTG(Era(BP)): Era", "MTG(Era(BP)): BP", "BP-Era", 
  "MTG(BP-Era): MTG", "MTG(BP-Era): BP-Era", "Clusters",
  "Clusters(BP-Era): Clusters", "Clusters(BP-Era): BP-Era"
)

write_csv(
  amova_table, file.path("results", "categorizations_AMOVA_variance.csv"), col_names = TRUE
)

