# import file paths and functions
source(file.path("src", "R", "file_paths.R"))
import::from(ape, "read.gff")
import::from(
  dplyr, "arrange", "do", "filter", "group_by", "left_join", "n", "rename",
  "select", "ungroup"
)
import::from(magrittr, "%>%")
import::from(parallel, "detectCores", "mclapply")
import::from(
  readr, "col_character", "col_double", "col_factor", "col_integer", "read_csv",
  "read_rds", "type_convert", "write_csv", "write_rds"
)
import::from(stringr, "str_c")
import::from(tibble, "as_tibble", "tibble")

chr_orders <- list(
  ABD = outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
    t() %>% as.vector(),
  ADB = outer(as.character(1:7), c("A", "D", "B"), paste, sep = "") %>%
    t() %>% as.vector(),
  BAD = outer(as.character(1:7), c("B", "A", "D"), paste, sep = "") %>%
    t() %>% as.vector(),
  BDA = outer(as.character(1:7), c("B", "D", "A"), paste, sep = "") %>%
    t() %>% as.vector(),
  DBA = outer(as.character(1:7), c("D", "B", "A"), paste, sep = "") %>%
    t() %>% as.vector(),
  DAB = outer(as.character(1:7), c("D", "A", "B"), paste, sep = "") %>%
    t() %>% as.vector()
)

# parse the GFF3 format file of the wheat 90K snp chip physical map positions
marker_aligns <- lapply(chr_orders[1], function (chr) {
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
    pos <- floor( (as.integer(feat[4]) + as.integer(feat[5]) ) / 2)
    return(c(name, ID, chrom, pos, coverage, per_id, use.names = F))
  }) %>% t()
}) %>% do.call(rbind, .)

# format the alignemnts, rename, filter, and split by marker
marker_aligns <- marker_aligns %>%
  as_tibble() %>%
  type_convert(
    col_type = list(
      col_character(), col_character(), col_character(), col_integer(),
      col_double(), col_double()
    )
  ) %>%
  rename(
    marker = V1, ID = V2, chrom = V3, pos = V4, coverage = V5, per_id = V6
  ) %>%
  filter(coverage >= 90, per_id >= 98) %>%
  split(.$marker)

# load the filtered gen map
marker_gen_pos <- read_rds(file.path(inter_markers, "pozniak_filtered_map.rds"))

# a function for finding the best alignment of each maker based on which
# linkage group it is assigned to
best_alignments <- function (marker_align) {
  gen_data <- marker_gen_pos[
    which(marker_gen_pos$marker == marker_align$marker[1]),
  ]
  if (nrow(gen_data)) {
    best_alignment <- marker_align[
      which(marker_align$chrom == gen_data$chrom),
    ]
    if (nrow(best_alignment) == 1) {
      return(best_alignment)
    }
  }
  return(tibble())
}

# remove at least one marker mapping to the same position in phys map
unique_marker <- function(markers) {
  if (nrow(markers) == 1) {
    return(markers[1, ])
  } else if (nrow(markers) == 2) {
    marker1 <- markers[1, 5:ncol(markers)]
    marker2 <- markers[2, 5:ncol(markers)]
    per_diff <- (sum(abs(marker1 - marker2), na.rm = TRUE) / 2) / 100
    if (per_diff > 0.01) {
        return(tibble())
    } else {
        return(markers[1, ])
    }
  } else {
    return(tibble())
  }
}

# format the genotype data into the proper format for snpgds format
genotypes <- read_csv(
  file.path(markers, "Jan_6_wheat_genotypes_curtis.csv")
) %>%
  select(-X1, -X3, -X4, -X5, -Name) %>%
  .[-1:-2, ] %>%
  rename(marker = X2) %>%
  replace(. == "C1", 0) %>%
  replace(. == "c1", 0) %>%
  replace(. == "C2", 2) %>%
  replace(. == "NC", 3)

blah <- lapply(chr_orders[1], function (chr_order) {
  levels(marker_gen_pos$chrom) <<- chr_order

  # identify the best alignments and combine with the gen map and genotypes
  markers_with_maps_genos <- mclapply(
      marker_aligns, best_alignments, mc.cores = detectCores()
    ) %>%
    do.call(rbind, .) %>%
    left_join(marker_gen_pos, by = c("marker", "chrom")) %>%
    rename(phys_pos = pos.x, gen_pos = pos.y) %>%
    select(marker, chrom, phys_pos, gen_pos) %>%
    group_by(chrom, phys_pos) %>%
    left_join(genotypes) %>%
    type_convert()

  # print the number of markers with identical phys pos
  markers_with_maps_genos %>%
    filter(n() >= 2) %>%
    nrow() %>%
    str_c("Num markers with identical postions: ", .) %>%
    print()

  # identify the number of the above markers we remove
  markers_with_maps_genos %>%
    filter(n() >= 2) %>%
    do(unique_marker(.)) %>%
    nrow() %>%
    str_c(
      "Of markers with identical postions num retained after pruning: ", .
    ) %>% print()

  # actually remove the markers
  markers_with_maps_genos <- markers_with_maps_genos %>%
    do(unique_marker(.)) %>%
    ungroup()
  nrow(markers_with_maps_genos) %>%
    str_c("Num markers remaining after pruning: ", .) %>%
    print()

  # write out the maps with genotypes
  write_rds(
    markers_with_maps_genos,
    file.path(
      inter_markers,
      str_c(
        chr_order[1], chr_order[2], chr_order[3], "maps_genotypes.rds",
        sep = "_"
      )
    )
  )
  write_csv(
    marker_gen_pos$marker[
      which(! marker_gen_pos$marker %in% markers_with_maps_genos$marker)
    ] %>% as.data.frame(),
    file.path(
      "results",
      str_c(
        chr_order[1], chr_order[2], chr_order[3], "filtered_markers.csv",
        sep = "_"
      )
    )
  )
})