# import file paths and functions
source("wheat-pgda/src/R/file_paths.R")
import::from(ape, "read.gff")
import::from(
  dplyr, "arrange", "do", "filter", "group_by", "left_join", "mutate", "n", "rename",
  "top_n", "select",  "ungroup"
)
library(ggplot2)
import::from(ggrepel, "geom_label_repel")
import::from(gridExtra, "grid.arrange")
import::from(magrittr, "%>%")
import::from(parallel, "detectCores", "mclapply")
import::from(pgda, "load_genes", "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(
  readr, "col_character", "col_double", "col_factor", "col_integer", "read_csv",
  "read_rds", "type_convert", "write_csv", "write_rds"
)
import::from(plyr, "rbind.fill")
import::from(tidyr, "drop_na")
import::from(stringr, "str_c")
import::from(tibble, "as_tibble", "tibble")

chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
    t() %>% as.vector()

phys_data <- snpgds_parse(phys_gds)

# parse the GFF3 format file of the wheat 90K snp chip physical map positions
marker_aligns1 <- lapply(chroms, function (chr) {
  chr_feats <- read.gff(
    file.path(
      markers, "90K_RefSeqv1_probe_alignments",
      str_c("Infinium90K-chr", chr, ".gff3")
    )
  )
  apply(chr_feats, 1, function (feat) {
    attributes <- strsplit(feat[9], split = ";")[[1]]
    name <- strsplit(attributes[2], split = "=")[[1]][2]
    ID <- strsplit(attributes[4], split = "=")[[1]][2]
    coverage <- strsplit(attributes[1], split = "=")[[1]][2]
    per_id <- strsplit(attributes[3], split = "=")[[1]][2]
    chrom <- substr(feat[1], 4, nchar(feat[1]))
    pos <- floor((as.integer(feat[4]) + as.integer(feat[5]) ) / 2)
    return(c(name, ID, chrom, pos, coverage, per_id, use.names = F))
  }) %>% t()
}) %>% do.call(rbind, .)

# format the alignemnts, rename, filter, and split by marker
marker_aligns2 <- marker_aligns1 %>%
  as_tibble() %>%
  type_convert(
    col_type = list(
      col_character(), col_character(), col_character(), col_integer(),
      col_double(), col_double()
    )
  ) %>%
  rename(
    marker = V1, ID = V2, chrom = V3, phys_pos = V4, coverage = V5, per_id = V6
  ) %>%
  filter(coverage >= 90, per_id >= 98)

# load the wang genetic map
wang_gen_map <- read_csv(file.path(markers, "wang_gen_map.csv"),
  na = "#N/A", col_names = c("marker", "chrom", "gen_pos"),
  col_types = list(
      col_character(), col_factor(NULL), col_double()
  )
) %>%
  arrange(chrom, gen_pos)

marker_aligns3 <- marker_aligns2 %>% left_join(wang_gen_map)
marker_aligns3 <- marker_aligns3[complete.cases(marker_aligns3), ]

max_phys_lengths <- span_by_chrom(
  marker_aligns3$chrom, marker_aligns3$phys_pos
) %>% max_lengths() / 1e6
max_gen_lengths <- span_by_chrom(
  marker_aligns3$chrom, marker_aligns3$gen_pos
) %>% max_lengths() / 100

landmarks <- rbind.fill(
  load_genes(
    file.path(blast, "selected_pheno.csv")
  ) %>% mutate(
    type = "Gene", base = 0.5, pos_mb = pos / 1e6
  ) %>%
    select(-pos),
  load_genes(
    file.path(blast, "selected_resi.csv")
  ) %>% mutate(
    type = "Gene", base = 0.5, pos_mb = pos / 1e6
  ) %>%
    select(-pos),
  cbind(
    id = "Centromere", type = "Centromere", base = 0.5,
    read_csv(file.path(intermediate, "centromeres.csv"))
  )
)

chrom_plots <- lapply(chroms, function (chrom) {
  chrom_landmarks <- landmarks[which(landmarks$chrom == chrom), ]

  rows <- which(marker_aligns3$chrom == chrom)

  kept_markers <- which(marker_aligns3$marker[rows] %in% phys_data$snp$id)

  colour_groups <- rep("Filtered", length(rows))
  colour_groups[kept_markers] <- "Kept"

  max_pos_cm <- max_gen_lengths[[
    ifelse(grepl("1", chrom), "one",
      ifelse(grepl("2", chrom), "two",
        ifelse(grepl("3", chrom), "three",
          ifelse(grepl("4", chrom), "four",
            ifelse(grepl("5", chrom), "five",
              ifelse(grepl("6", chrom), "six", "seven")
            )
          )
        )
      )
    )
  ]]

  max_pos_mb <- max_phys_lengths[[
    ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
  ]]

  ggplot() +
    expand_limits(x = c(0, max_pos_mb), y = c(0, max_pos_cm)) +
    geom_point(
      aes(
        marker_aligns3[rows, ]$phys_pos, marker_aligns3[rows, ]$gen_pos,
        colour = colour_groups
      ),
      size = 2
    ) +
    labs(x = "Position in Mb", y = "Position in cM", title = chrom) +
    geom_vline(
      aes(
        xintercept = chrom_landmarks$pos_mb,
        linetype = chrom_landmarks$type
      ), size = 1
    ) +
    geom_label_repel(
      aes(
        chrom_landmarks$pos_mb, max_pos_cm,
        label = chrom_landmarks$id
      ), size = 2, seed = 101, inherit.aes = FALSE
    ) +
    scale_linetype("Landmarks")
})  


# plot
png(
  file.path("results", "phys_pos_vs_gen_pos_no_filtering.png"),
  family = "Times New Roman", width = 1800, height = 2800, pointsize = 5,
  units = "mm", res = 96
)
grid.arrange(
  grobs = chrom_plots, nrow = 7, ncol = 3
)
dev.off()