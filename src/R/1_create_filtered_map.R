# load import needed functions and file paths
source("wheat-pgda/src/R/file_paths.R")
import::from(dplyr, "anti_join", "arrange", "bind_rows", "semi_join")
import::from(magrittr, "%>%")
import::from(
  readr, "col_character", "col_double", "col_factor", "read_csv", "write_rds"
)

# load the pozniak genetic map
poz_gen_map <- read_csv(file.path(markers, "pozniak_gen_map.csv"),
  na = "#N/A", col_names = c("marker", "chrom", "pos"),
  col_types = list(
    col_character(), col_factor(NULL), col_double()
  )
) %>%
  na.omit() %>%
  arrange(chrom, pos)

# load the wang genetic map
wang_gen_map <- read_csv(file.path(markers, "wang_gen_map.csv"),
  na = "#N/A", col_names = c("marker", "chrom", "pos"),
  col_types = list(
      col_character(), col_factor(NULL), col_double()
  )
) %>%
  arrange(chrom, pos)

# filter the genetic maps so only markers present in the poz gen map are kept
# with those on a different chrom in the wang genetic map removed
poz_semi <- semi_join(poz_gen_map, wang_gen_map, by = c("marker", "chrom"))
poz_anti <- anti_join(poz_gen_map, wang_gen_map, by = c("marker", "chrom")) %>%
  anti_join(wang_gen_map, by = "marker")
poz_filtered <- poz_semi %>%
  bind_rows(poz_anti) %>%
  arrange(marker, chrom, pos)

# write the filtered genetic map out
write_rds(poz_filtered, file.path(inter_markers, "pozniak_filtered_map.rds"))

