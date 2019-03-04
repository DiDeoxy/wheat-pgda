# load import needed functions and file paths
source(file.path("src", "R", "file_paths.R"))
import::from(dplyr, "anti_join", "arrange", "bind_rows", "semi_join")
import::from(magrittr, "%>%")
import::from(
  readr, "col_character", "col_double", "col_factor", "read_csv", "write_rds"
)

# load the pozniak genetic map
poz_gen_map <- read_csv(file.path(markers, "pozniak_gen_map.csv"),
  na = "#N/A", col_names = c("marker", "group", "pos"),
  col_types = list(
    col_character(), col_factor(NULL), col_double()
  )
) %>%
  na.omit() %>%
  arrange(group, pos)

# load the wang genetic map
wang_gen_map <- read_csv(file.path(markers, "wang_gen_map.csv"),
  na = "#N/A", col_names = c("marker", "group", "pos"),
  col_types = list(
      col_character(), col_factor(NULL), col_double()
  )
) %>%
  arrange(group, pos)

# filter the genetic maps so only markers present in the poz gen map are kept
# with those having a different position in the wang genetic map removed
poz_semi <- semi_join(poz_gen_map, wang_gen_map, by = c("marker", "group"))
poz_anti <- anti_join(poz_gen_map, wang_gen_map, by = c("marker", "group")) %>%
  anti_join(wang_gen_map, by = "marker")
poz_filtered <- poz_semi %>%
  bind_rows(poz_anti) %>%
  arrange(marker, group, pos)

# write the filtered genetic map out
write_rds(poz_filtered, file.path(inter_markers, "pozniak_filtered_map.rds"))

