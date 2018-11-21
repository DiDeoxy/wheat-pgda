library(tidyverse)

# load the pozniak genetic map
poz_gen_map <- read_csv("Data/Raw/Maps/pozniak_gen_map.csv",
    na = "#N/A",
    col_names = c("marker", "group", "pos"),
    col_types = list(
        col_character(), col_factor(NULL), col_double()
    )
) %>%
    na.omit() %>%
    arrange(group, pos)

# load the wang genetic map
wang_gen_map <- read_csv("Data/Raw/Maps/wang_gen_map.csv",
    na = "#N/A",
    col_names = c("marker", "group", "pos"),
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
write_rds(
    poz_filtered,
    "Data/Intermediate/Maps/pozniak_filtered_map.rds"
)