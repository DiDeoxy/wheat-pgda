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
  "read_rds", "type_convert", "write_rds"
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
marker_aligns <- lapply(chr_orders$ABD, function (chr) {
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

# a function of rfinding the best alignment of each maker based on which
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

# load the filtered gen map
marker_gen_pos <- read_rds(file.path(inter_markers, "pozniak_filtered_map.rds"))

# add the genes positons to the regions table
snp_data <- wheat_data$snp %>%
  mutate(pos_mb = pos / 1e6) %>%
  arrange(chrom, pos_mb) %>%
  select(-pos) %>%
  add_column(josts_d := read_rds(josts_d)) %>%
  cbind(mjafs_by_pop)

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

snp_data <- snp_data %>%
  rowwise() %>%
  mutate(class =
    ifelse(
      josts_d > top_quartile,
      c("CHRSD", "CHRWD", "CSWSD")[
        which.max(
          c(
            sum(abs(chrs - chrw), abs(chrs - csws)),
            sum(abs(chrw - chrs), abs(chrw - csws)),
            sum(abs(csws - chrs), abs(csws - chrw))
          )
        )
      ],
      "None"
    )
  )

# load the gene positions
all_genes <- rbind.fill(
  load_genes(
    file.path(blast, "selected_pheno.csv"), base = 1
  ) %>% mutate(gene_type = "Phenotype Genes", pos_mb = pos / 1e6) %>%
    select(-pos),
  load_genes(
    file.path(blast, "selected_resi.csv"), base = 1
  ) %>% mutate(gene_type = "Resistance Genes", pos_mb = pos / 1e6) %>%
    select(-pos)
)

snp_data <- snp_data %>%
  rbind.fill(all_genes) %>%
  mutate(
    type = pmin(gene_type, class, na.rm = TRUE),
    jdb = pmin(josts_d, base, na.rm = TRUE)
  ) %>%
  select(-c(base, josts_d, gene_type, class)) %>%
  spread(type, jdb) %>%
  gather(value_type, values, -c(id, chrom, pos_mb)) %>%
  arrange(chrom, pos_mb)

  # print the number of markers with identical phys pos
  markers_with_maps_genos %>%
    filter(n() >= 2) %>%
    nrow() %>%
    str_c("Num markers with identical postions: ", .) %>%
    print()

genes_nearby_markers <- lapply(1:nrow(gene_ranges), function (row) {
  add_column(
    snp_data[
      which(
        snp_data$chrom == gene_ranges[row, ]$chrom &
        snp_data$pos_mb >= gene_ranges[row, ]$window_start &
        snp_data$pos_mb <= gene_ranges[row, ]$window_end
      ),
    ],
    group = gene_ranges[row, ]$genes,
    window_start = gene_ranges[row, ]$window_start,
    window_end = gene_ranges[row, ]$window_end
  )
})

genes_nearby_markers[[40]]

legend_title <- "MJAF, Jost's D,\nand Gene Type"
lables <- c(
  "CHRS MJAF", "CHRSD Jost's D", "CHRW MJAF", "CHRW Jost's D", "CSWS MJAF",
  "CSWSD Jost's D", "None Jost's D", "Phenotype Genes", "Resistance Genes"
)
lapply(genes_nearby_markers, function(gene_data) {
  if (nrow(gene_data)) {
    png(
      file.path(
        zoomed_marker_plots,
        str_c(
          str_c(
            gene_data$chrom[1],
            gene_data$window_start[1], gene_data$window_end[1],
            gene_data$group[1],
            sep = "_"
          ),
          ".png"
        )
      ),
      family = "Times New Roman",
      width = 320, height = 240,
      units = "mm", res = 192
    )
    print(gene_data %>%
      ggplot(aes(pos_mb, values)) +
      ylim(0, 1) +
      xlim(gene_data$window_start[1], gene_data$window_end[1]) +
      geom_point(aes(colour = value_type, shape = value_type), size = 2) +
      geom_hline(yintercept = top_quartile) +
      geom_text_repel(
        aes(
          label = ifelse(
            value_type == "Phenotype Genes" | value_type == "Resistance Genes",
            id, ""
          )
        )
      ) +
      scale_colour_manual(
        legend_title, labels = lables,
        values = colour_set[c(1, 1, 2, 2, 4, 4, 20, 15, 19)]
      ) +
      scale_shape_manual(
        legend_title, labels = lables,
        values = c(16, 11, 16, 11, 16, 11, 11, 25, 25)
      ) +
      labs(
        x = "Position in Mb", y = "MJAF and Jost's D",
        title = str_c(gene_data$chrom[1], gene_data$group[1], sep = "_")
      )
    )
    dev.off()
  }
})