library(DEMEtics)

for (gene in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")[1]) {
  blah <- Gst.Nei(
    str_c("Data\\Intermediate\\DEMEtics\\", gene, "_genind.tsv"),
  format.table = FALSE, pm = "overall", statistics = "p", bt = 100)
}