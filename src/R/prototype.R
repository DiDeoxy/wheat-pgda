source(file.path("src", "R", "file_paths.R"))
library(pgda)
library(magrittr)

x_data <- list(
  phys_data = snpgds_parse(phys_gds), gen_data = snpgds_parse(gen_gds),
  ld_phys_data = snpgds_parse(ld_phys_gds),
  ld_gen_data = snpgds_parse(ld_gen_gds)
)

for (name in names(x_data)) {
  print(name)

  total_span <- span_by_chrom(
    x_data[[name]]$snp$chrom, x_data[[name]]$snp$pos, diff = TRUE
  ) %>% sum()
  print(total_span)
  print(total_span / nrow(x_data[[name]]$snp))
  print(
    coverage_by_chrom(x_data[[name]]$snp$chrom, x_data[[name]]$snp$pos) %>%
    sum() / total_span
  )
}