setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading and ordering the SNP data genetic map
pozniak_gen_map <- read.csv("Data\\Raw\\Maps\\pozniak_gen_map.csv",
                        header = F, col.names = c("Name", "Contig", "Position"), row.names = 1, stringsAsFactors = F)
pozniak_gen_map <- pozniak_gen_map[-which(pozniak_gen_map[,1] == "#N/A"),]
pozniak_gen_map <- pozniak_gen_map[order(pozniak_gen_map$Contig, as.numeric(pozniak_gen_map$Position)),]
dim(pozniak_gen_map)
## loading and ordering the S13 genetic map
wang_gen_map <- read.csv("Data\\Raw\\Maps\\wang_gen_map.csv",
                        header = F, col.names = c("Name", "Contig", "Position"), row.names = 1, stringsAsFactors = F)
wang_gen_map <- wang_gen_map[order(wang_gen_map$Contig, as.numeric(wang_gen_map$Position)),]
dim(wang_gen_map)
## finding the overlap by name between the two maps
pozniak_overlap <- pozniak_gen_map[row.names(pozniak_gen_map) %in% row.names(wang_gen_map),]
wang_overlap <- wang_gen_map[row.names(wang_gen_map) %in% row.names(pozniak_gen_map),]
## comparing the order of markers between the two
overlap_order <- cbind(row.names(pozniak_overlap), row.names(wang_overlap))

## identifying the overlap markers that are also on the same chromosome
agree_probes <- list()
for (name in row.names(wang_overlap)){
  if (wang_gen_map[name,1] == pozniak_gen_map[name,1]) {
    agree_probes[[name]] <- c(name, pozniak_gen_map[name,1], pozniak_gen_map[name,2])
  }
}
agree_probes <- do.call(rbind.data.frame, agree_probes)
names(agree_probes) <- c("Name", "Contig", "Position")
row.names(agree_probes) <- agree_probes[,1]
agree_probes[,1] <- NULL
dim(agree_probes)
## identifying the overlap markers that are on different chromosomes
disagree_probes <- list()
for (name in row.names(wang_overlap)){
  if (wang_gen_map[name,1] != pozniak_gen_map[name,1]) {
    disagree_probes[[name]] <- c(name, wang_gen_map[name,1], pozniak_gen_map[name,1])
  }
}
disagree_probes <- do.call(rbind.data.frame, disagree_probes)
names(disagree_probes) <- c("Name", "Contig_wang", "Position_pozniak")
row.names(disagree_probes) <- disagree_probes[,1]
disagree_probes[,1] <- NULL
dim(disagree_probes)
## identifying the unique sets of markers in each 
pozniak_unique <- pozniak_gen_map[! row.names(pozniak_gen_map) %in% row.names(wang_gen_map),]
wang_unique <- wang_gen_map[! row.names(wang_gen_map) %in% row.names(pozniak_gen_map),]

## Different maps
pozniak_filtered_map <- rbind(agree_probes, pozniak_unique, stringsAsFactors = F)
pozniak_filtered_map <- pozniak_filtered_map[order(pozniak_filtered_map$Contig, as.numeric(pozniak_filtered_map$Position)),]
dim(pozniak_filtered_map)
wang_filtered_map <- rbind(agree_probes, wang_unique, stringsAsFactors = F)
wang_filtered_map <- wang_filtered_map[order(wang_filtered_map$Contig, as.numeric(wang_filtered_map$Position)),]

combined_filtered_map <- rbind(agree_probes, wang_unique, pozniak_unique, stringsAsFactors = F)
combined_filtered_map <- combined_filtered_map[order(combined_filtered_map$Contig, as.numeric(combined_filtered_map$Position)),]

## writing the maps out to file
write.csv(pozniak_filtered_map, "Data\\Intermediate\\Maps\\Phys\\pozniak_filtered_map.csv")
write.csv(wang_filtered_map, "Data\\Intermediate\\Maps\\Phys\\wang_filtered_map.csv")
write.csv(combined_filtered_map, "Data\\Intermediate\\Maps\\Phys\\combined_filtered_map.csv")
write.csv(agree_probes, "Data\\Intermediate\\Maps\\Phys\\overlap_map.csv")
