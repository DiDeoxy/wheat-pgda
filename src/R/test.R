source(file.path("src", "R", "file_paths.R"))
import::from(pgda, "snpgds_parse")
import::from(magrittr, "%>%")

wheat_data <- snpgds_parse(phys_gds)

# table(paste0(wheat_data$sample$annot$bp, wheat_data$sample$annot$era)) %>% sort()
# table(paste0(wheat_data$sample$annot$bp, wheat_data$sample$annot$era)) %>% length()

# table(data.frame(wheat_data$sample$annot$bp, wheat_data$sample$annot$era))

# table(wheat_data$sample$annot$pheno)
# table(wheat_data$sample$annot$mc)
# table(wheat_data$sample$annot$pheno, wheat_data$sample$annot$era)
# table(wheat_data$sample$annot$mc, wheat_data$sample$annot$era)

