source(file.path("src", "R", "file_paths.R"))
import::from(pgda, "coverage_by_chrom", "snpgds_parse")
import::from(magrittr, "%>%")
import::from(stringr, "str_c", "str_wrap")
import::from(tibble, "tibble")

wheat_data <- snpgds_parse(phys_gds)
(coverage_by_chrom(wheat_data$snp$chrom, wheat_data$snp$pos) %>% sum()) / (wheat_data$chrom_lengths %>% sum())
wheat_data$chrom_lengths[seq(3, 21, 3)] %>% sum() / wheat_data$chrom_lengths %>% sum()
chrom_num_markers <- by(wheat_data$snp, wheat_data$snp$chrom, nrow) %>% as.list() %>% unlist()
chrom_num_markers[seq(3, 21, 3)] %>% sum() / chrom_num_markers %>% sum()

wheat_data <- snpgds_parse(gen_gds)

(coverage_by_chrom(wheat_data$snp$chrom, wheat_data$snp$pos) %>% sum()) / (wheat_data$chrom_lengths %>% sum())