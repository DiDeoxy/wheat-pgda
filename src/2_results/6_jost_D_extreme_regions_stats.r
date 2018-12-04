library(adegenet)
library(tidyverse)

source("src/R_functions/funcs_gds_parse_create.R")

wheat_data <- parse_gds("phys_subset_sample")

################################################################################
# multi

# 1B
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 2 &
      wheat_data$snp$pos_mb >= 678.635327 &
      wheat_data$snp$pos_mb <= 680.17122
  ),
] %>% nrow()

# 3A 1
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 7 &
      wheat_data$snp$pos_mb >= 100.881838 &
      wheat_data$snp$pos_mb <= 106.435192
  ),
] %>% nrow()

# 3A 2
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 7 &
      wheat_data$snp$pos_mb >= 540.669124 &
      wheat_data$snp$pos_mb <= 543.712443
  ),
] %>% nrow()

# 3D
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 7 &
      wheat_data$snp$pos_mb >= 570.80153 &
      wheat_data$snp$pos_mb <= 571.550497
  ),
] %>% nrow()


################################################################################
comp_genind <- read_rds("Data/Intermediate/Adegenet/chrs_chrw_genind.rds")

# 2D
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 6 &
      wheat_data$snp$pos_mb >= 76.007955 &
      wheat_data$snp$pos_mb <= 78.793706
  ),
] %>% nrow()

markers <- c(least = "Kukri_c16557_855", most = "BS00082605_51")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 3B
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 8 &
      wheat_data$snp$pos_mb >= 488.187657 &
      wheat_data$snp$pos_mb <= 507.510674
  ),
] %>% nrow()

markers <- c(least = "Excalibur_c46052_695", most = "BS00030430_51")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 7B
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 20 &
      wheat_data$snp$pos_mb >= 511.064039 &
      wheat_data$snp$pos_mb <= 517.327995
  ),
] %>% nrow()

markers <- c(least = "Tdurum_contig53901_177", most = "Excalibur_c4556_113")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

################################################################################
comp_genind <- read_rds("Data/Intermediate/Adegenet/chrs_csws_genind.rds")

# 1D
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 3 &
      wheat_data$snp$pos_mb >= 407.875138 &
      wheat_data$snp$pos_mb <= 416.456286
  ),
] %>% nrow()

markers <- c(least = "BS00030252_51", most = "Ra_c21676_178")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 2A
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 4 &
      wheat_data$snp$pos_mb >= 748.968009 &
      wheat_data$snp$pos_mb <= 752.092859
  ),
] %>% nrow()

markers <- c(both = "CAP7_c4056_108")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 3A 1
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 7 &
      wheat_data$snp$pos_mb >= 519.273991 &
      wheat_data$snp$pos_mb <= 525.200911
  ),
] %>% nrow()

markers <- c(least = "BS00026189_51", most = "Ku_c56370_1155")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 3A 2
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 7 &
      wheat_data$snp$pos_mb >= 625.239797 &
      wheat_data$snp$pos_mb <= 627.825027
  ),
] %>% nrow()

markers <- c(least = "wsnp_Ex_c9483_15722127", most = "wsnp_Ex_c8884_14841846")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 3B
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 8 &
      wheat_data$snp$pos_mb >= 682.905538 &
      wheat_data$snp$pos_mb <= 700.793524
  ),
] %>% nrow()

markers <- c(least = "tplb0062a17_679", most = "wsnp_CAP11_c2309_1201554")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 5B
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 14 &
      wheat_data$snp$pos_mb >= 571.475214 &
      wheat_data$snp$pos_mb <= 573.805418
  ),
] %>% nrow()

markers <- c(both = "RAC875_c30011_426")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 7A
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 19 &
      wheat_data$snp$pos_mb >= 225.356544 &
      wheat_data$snp$pos_mb <= 244.470445
  ),
] %>% nrow()

markers <- c(both = "wsnp_JD_c18167_16742264", most = "RAC875_rep_c97817_197")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 7B
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 20 &
      wheat_data$snp$pos_mb >= 739.38387 &
      wheat_data$snp$pos_mb <= 742.591455
  ),
] %>% nrow()

markers <- c(least = "Tdurum_contig53901_177", most = "Excalibur_c4556_113")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

################################################################################
comp_genind <- read_rds("Data/Intermediate/Adegenet/csws_chrw_genind.rds")

# 1D
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 3 &
      wheat_data$snp$pos_mb >= 11.455911 &
      wheat_data$snp$pos_mb <= 12.53452
  ),
] %>% nrow()

markers <- c(both = "BS00087783_51")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 2B 1
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 5 &
      wheat_data$snp$pos_mb >= 72.575726 &
      wheat_data$snp$pos_mb <= 77.291044
  ),
] %>% nrow()

markers <- c(least = "Kukri_rep_c101093_572", most = "RAC875_c24914_462")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 2B 2
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 5 &
      wheat_data$snp$pos_mb >= 140.750258 &
      wheat_data$snp$pos_mb <= 192.364853
  ),
] %>% nrow()

markers <- c(least = "Kukri_rep_c102367_278", most = "Excalibur_c23185_155")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 3B
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 8 &
      wheat_data$snp$pos_mb >= 116.147522 &
      wheat_data$snp$pos_mb <= 133.201248
  ),
] %>% nrow()

markers <- c(both = "Kukri_c39386_282")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 4B
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 11 &
      wheat_data$snp$pos_mb >= 652.253149 &
      wheat_data$snp$pos_mb <= 657.526033
  ),
] %>% nrow()

markers <- c(both = "BS00022055_51")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

# 6B
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 17 &
      wheat_data$snp$pos_mb >= 199.718577 &
      wheat_data$snp$pos_mb <= 276.224891
  ),
] %>% nrow()

markers <- c(least = "BS00025302_51", most = "RFL_Contig1492_342")

for (marker in markers) {
  print(marker)
  tibble(
    pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}



# num extreme markers in significant loci
# CHRS
59 + 12 + 11 + 20 + 14 + 25 + 5 + 11 + 2 + 28 + 19 + 9
# CHRW
65 + 12 + 20 + 14 + 5 + 22 + 54 + 19
# CSWS
59 + 65 + 11 + 25 + 5 + 5 + 22 + 11 + 2 + 54 + 28 + 9
