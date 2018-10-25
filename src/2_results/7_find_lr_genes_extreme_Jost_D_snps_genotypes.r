library(adegenet)

# find the underlying genetoypes of the extreme markers near genes

# Lr10
comp_genind <- read_rds(str_c("Data\\Intermediate\\Adegenet\\Lr10_genind.rds"))

markers <- c(
    "tplb0062i18_102"
)

for (marker in markers) {
    print(marker)
    cbind(
        allele = comp_genind$tab[, str_c(marker, ".A")], comp_genind$strata
    ) %>% table() %>% print()
}

# Lr21
comp_genind <- read_rds(str_c("Data\\Intermediate\\Adegenet\\Lr21_genind.rds"))

markers <- c("wsnp_Ex_c1358_2600929", "wsnp_Ex_c1358_2602235")

for (marker in markers) {
    print(marker)
    cbind(
        allele = comp_genind$tab[, str_c(marker, ".A")], comp_genind$strata
    ) %>% table() %>% print()
}

# lr1
comp_genind <- read_rds(str_c("Data\\Intermediate\\Adegenet\\Lr1_genind.rds"))

markers <- c(
    "tplb0031c19_721", "RAC875_c8100_245", "D_GB5Y7FA01D4CBK_54",
    "RAC875_c34968_514", "Kukri_c21384_1475", "BS00012069_51",
    "Excalibur_c24145_1643", "wsnp_Ex_c24145_33394561"
)

for (marker in markers) {
    print(marker)
    cbind(
        allele = comp_genind$tab[, str_c(marker, ".A")], comp_genind$strata
    ) %>% table() %>% print()
}

# Lr34
comp_genind <- read_rds(str_c("Data\\Intermediate\\Adegenet\\Lr34_genind.rds"))

markers <- c("Kukri_c92151_216", "Kukri_c32845_116", "TA002473-0717")

for (marker in markers) {
    print(marker)
    cbind(
        allele = comp_genind$tab[, str_c(marker, ".A")], comp_genind$strata
    ) %>% table() %>% print()
}