library(tidyverse)
library(magrittr)
library(SNPRelate)
library(plyr)

gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")

makeGenesContigsDataframe <- function (file_in, file_out) {
    genes <- read_csv(file_in, 
                      col_names = c("name", "ID", "contig", "pos", "sleng", 
                                    "aleng", "%id")) %>%
             select(name, contig, pos) %>%
             arrange(contig, pos)
    
    chroms <-as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste,
                               sep="")))
    chroms_part2 <- paste0(chroms, "_part2")
    part2_start <- tibble(contig = chroms_part2,
                          start = c(471304005, 438720154, 452179604, 462376173, 
                                    453218924, 462216879, 454103970, 448155269,
                                    476235359, 452555092, 451014251, 451004620,
                                    453230519, 451372872, 451901030, 452440856, 
                                    452077197, 450509124, 450046986, 453822637,
                                    453812268))
    genes_crctd <- genes %>% 
    left_join(part2_start) %<>% 
    mutate(pos = rowSums(cbind(.$pos, .$start), na.rm = TRUE)) %<>%
    mutate(contig = substr(.$contig, 1, 2)) %>%
    select(name, contig, pos)

    genes_contigs <- tibble()
    for (i in unique(snp_chrom)) {
        chrom_genes <- genes_crctd %>% filter(contig == chroms[i])
        remaining <- length(which(snp_chrom == i)) - nrow(chrom_genes)
        empty <- tibble(name = rep(NA, remaining),
                        contig = rep(NA, remaining),
                        pos = rep(NA, remaining))
        genes_contigs <- genes_contigs %>% rbind(chrom_genes, empty)
    }
    
    write_rds(genes_contigs, path = file_out)
}

makeGenesContigsDataframe(
    "Data\\Intermediate\\Aligned_genes\\top_main_genes_selected.csv", 
    "Data\\Intermediate\\Aligned_genes\\top_main_genes_selected_contigs.rds")
makeGenesContigsDataframe(
    "Data\\Intermediate\\Aligned_genes\\top_resistance_genes_selected.csv",
    paste0("Data\\Intermediate\\Aligned_genes\\top_resistance_genes_selected",
    "_contigs.rds"))
makeGenesContigsDataframe(
    "Data\\Intermediate\\Aligned_genes\\known_genes_groups.csv",
    "Data\\Intermediate\\Aligned_genes\\known_genes_groups.rds")
