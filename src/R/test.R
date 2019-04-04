source(file.path("src", "R", "file_paths.R"))
import::from(pgda, "snpgds_parse")

wheat_data <- snpgds_parse(gen_gds)

length(wheat_data$snp$pos) / length(unique(wheat_data$snp$pos))