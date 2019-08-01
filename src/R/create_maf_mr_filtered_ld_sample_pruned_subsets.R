# import file paths and functions
source(file.path("src", "R", "file_paths.R"))
import::from(dplyr, "as_tibble")
import::from(igraph, "graph_from_edgelist", "max_cliques")
import::from(magrittr, "%>%")
import::from(pgda, "snpgds_parse", "snpgds_sample_subset", "snpgds_snp_subset")
import::from(
  SNPRelate, "snpgdsClose", "snpgdsIBS", "snpgdsLDpruning", "snpgdsOpen",
  "snpgdsSelectSNP"
)
import::from(stringr, "str_c")

# create ld pruned set of markes
wheat_gds <- snpgdsOpen(file.path(gds, "phys.gds"))
set.seed(1000)
kept_id <- snpgdsLDpruning(
  wheat_gds, autosome.only = FALSE, maf = 0.05, missing.rate = 0.10, 
  slide.max.bp = 1e7, ld.threshold = 0.7
) %>% unlist()
snpgdsClose(wheat_gds)

wheat_data <- snpgds_parse(file.path(gds, "phys.gds"))
kept_index <- match(kept_id, wheat_data$snp$id) %>% sort()
snpgds_snp_subset(
  wheat_data, file.path(gds, "maf_mr_filtered_ld_pruned_phys.gds"), kept_index
)

wheat_data <- snpgds_parse(file.path(gds, "gen.gds"))
kept_index <- match(kept_id, wheat_data$snp$id) %>%sort()
snpgds_snp_subset(
  wheat_data, file.path(gds, "maf_mr_filtered_ld_pruned_gen.gds"), kept_index
)

# eliminate those individuals that show identity by state
# (IBS, fractional identity) greater than 0.99
wheat_gds <- snpgdsOpen(file.path(gds, "maf_mr_filtered_ld_pruned_phys.gds"))
IBS <- snpgdsIBS(wheat_gds, autosome.only = FALSE)
snpgdsClose(wheat_gds)

pairs <- which(IBS$ibs >= 0.99, arr.ind = TRUE)
pairs <- cbind(IBS$sample.id[pairs[, 1]], IBS$sample.id[pairs[, 2]])

indices <- vector()
for (i in 1:dim(pairs)[1]) {
  if (pairs[i, 1] == pairs[i, 2]) {
    indices <- c(indices, i)
  }
}
pairs <- pairs[-indices, ]

# print out a table of the maximal cliques constructed from the pairs
graph_from_edgelist(pairs) %>%
  max_cliques() %>%
  lapply(names) %>%
  lapply(`length<-`, max(lengths(.))) %>%
  do.call(rbind, .) %>%
  cbind(str_c("Clique ", 1:nrow(.)), .) %>%
  as_tibble() %>%
  write.table(
    file.path("results", "clique_table.csv"), sep = ",",
    row.names = FALSE, quote = FALSE,
    col.names = c("Clique", str_c("Cultivar ", 1:(ncol(.) - 1)))
  )

# used clique_table.csv to identify cultivars to prune from each clique
NILs <- c(
  "CDC Stanley 1", "PT754", "SWS349", "Somerset 1","Stettler 1",
  "PT434", "BW811", "AC Minto 1", "Avocet 1", "BW275 1",
  "BW395", "PT616", "BW427 1", "BW492", "BW948",
  "Carberry 1", "AC Reed 1", "SWS87", "SWS345", "SWS241",
  "SWS410", "SWS390", "SWS408"
)

# mr pruned phys map
wheat_data <- snpgds_parse(file.path(gds, "phys.gds"))
sample_index <- match(NILs, wheat_data$sample$id)
# create gds object without the NILs
snpgds_sample_subset(wheat_data, phys_gds, sample_index)

# mr pruned gen map
wheat_data <- snpgds_parse(file.path(gds, "gen.gds"))
sample_index <- match(NILs, wheat_data$sample$id)
snpgds_sample_subset(wheat_data, gen_gds, sample_index)

# mr pruned phys map
wheat_data <- snpgds_parse(file.path(gds, "maf_mr_filtered_ld_pruned_phys.gds"))
sample_index <- match(NILs, wheat_data$sample$id)
# create gds object without the NILs
snpgds_sample_subset(wheat_data, ld_phys_gds, sample_index)

# mr pruned phys map
wheat_data <- snpgds_parse(file.path(gds, "maf_mr_filtered_ld_pruned_gen.gds"))
sample_index <- match(NILs, wheat_data$sample$id)
# create gds object without the NILs
snpgds_sample_subset(wheat_data, ld_gen_gds, sample_index)