library(poppr)

## Breeding program 
# Per Locus data re-formatiing
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_bp.RData")
wheat.poppr.seploc <- seploc(wheat.poppr.bp)
for (i in 1:length(wheat.poppr.seploc)){
  strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.bp)
}
# per locus amova
wheat.amova <- lapply(wheat.poppr.seploc, poppr.amova, hier = ~bp, cutoff = 0.11)
# saving the amova object
save(wheat.amova, file = "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\AMOVA\\amova_bp.RData")


## Origin
# Per Locus data re-formatting
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_bp.RData")
wheat.poppr.seploc <- seploc(wheat.poppr.or)
for (i in 1:length(wheat.poppr.seploc)){
  strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.or)
}
# per locus amova
wheat.amova <- lapply(wheat.poppr.seploc, poppr.amova, hier = ~origin, cutoff = 0.11)
# saving the amova object
save(wheat.amova, file = "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\AMOVA\\amova_or.RData")

# ## strength
# # Per Locus data re-formatiing
# load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_st.RData")
# wheat.poppr.seploc <- seploc(wheat.poppr.st)
# for (i in 1:length(wheat.poppr.seploc)){
#   strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.st)
# }
# # per locus amova
# wheat.amova <- lapply(wheat.poppr.seploc, poppr.amova, hier = ~strength, cutoff = 0.11)
# # saving the amova object
# save(wheat.amova, file = "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\AMOVA\\amova_st.RData")
# 
# 
# ## colour
# # Per Locus data re-formatiing
# load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_co.RData")
# wheat.poppr.seploc <- seploc(wheat.poppr.co)
# for (i in 1:length(wheat.poppr.seploc)){
#   strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.co)
# }
# # per locus amova
# wheat.amova <- lapply(wheat.poppr.seploc, poppr.amova, hier = ~colour, cutoff = 0.11)
# # saving the amova object
# save(wheat.amova, file = "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\AMOVA\\amova_co.RData")

## srength + colour
# Per Locus data re-formatiing
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_hr_sw.RData")
wheat.poppr.seploc <- seploc(wheat.poppr.hr.sw)
for (i in 1:length(wheat.poppr.seploc)){
  strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.hr.sw)
}
# per locus amova
wheat.amova <- lapply(wheat.poppr.seploc, poppr.amova, hier = ~hr.sw, cutoff = 0.11)
# saving the amova object
save(wheat.amova, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\AMOVA\\amova_hr_sw.RData")

## season
# Per Locus data re-formatiing
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_se.RData")
wheat.poppr.seploc <- seploc(wheat.poppr.se)
for (i in 1:length(wheat.poppr.seploc)){
  strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.se)
}
# per locus amova
wheat.amova <- lapply(wheat.poppr.seploc, poppr.amova, hier = ~season, cutoff = 0.11)
# saving the amova object
save(wheat.amova, file = "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\AMOVA\\amova_se.RData")


## designation
# Per Locus data re-formatiing
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_de.RData")
wheat.poppr.seploc <- seploc(wheat.poppr.de)
for (i in 1:length(wheat.poppr.seploc)){
  strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.de)
}
# per locus amova
wheat.amova <- lapply(wheat.poppr.seploc, poppr.amova, hier = ~designation, cutoff = 0.11)
# saving the amova object
save(wheat.amova, file = "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\AMOVA\\amova_de.RData")


## market class
# Per Locus data re-formatiing
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_mc.RData")
wheat.poppr.seploc <- seploc(wheat.poppr.mc)
for (i in 1:length(wheat.poppr.seploc)){
  strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.mc)
}
# per locus amova
wheat.amova <- lapply(wheat.poppr.seploc, poppr.amova, hier = ~mc, cutoff = 0.11)
# saving the amova object
save(wheat.amova, file = "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\AMOVA\\amova_mc.RData")

## pamk
# Per Locus data re-formatiing
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_pk.RData")
wheat.poppr.seploc <- seploc(wheat.poppr.pk)
for (i in 1:length(wheat.poppr.seploc)){
  strata(wheat.poppr.seploc[[i]]) <- strata(wheat.poppr.pk)
}
# per locus amova
wheat.amova <- lapply(wheat.poppr.seploc, poppr.amova, hier = ~pamk, cutoff = 0.11)
# saving the amova object
save(wheat.amova, file = "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\AMOVA\\amova_pk.RData")