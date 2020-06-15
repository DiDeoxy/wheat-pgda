# # identify those snps within the extended haplotypes
# snp_index_1A <- which(wheat_data_phys$snp$chrom == "1A")
# haplo_id_1A <- wheat_data_phys$snp$id[snp_index_1A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_1A] > 70 &
#     wheat_data_phys$snp$pos_mb[snp_index_1A] < 300)
# ]
# snp_index_2A <- which(wheat_data_phys$snp$chrom == "2A")
# haplo_id_2A <- wheat_data_phys$snp$id[snp_index_2A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_2A] > 210 &
#     wheat_data_phys$snp$pos_mb[snp_index_2A] < 470)
# ]
# snp_index_4A <- which(wheat_data_phys$snp$chrom == "4A")
# haplo_id_4A <- wheat_data_phys$snp$id[snp_index_4A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_4A] > 230 &
#     wheat_data_phys$snp$pos_mb[snp_index_4A] < 460)
# ]
# snp_index_5B <- which(wheat_data_phys$snp$chrom == "5B")
# haplo_id_5B <- wheat_data_phys$snp$id[snp_index_5B][
#   which(wheat_data_phys$snp$pos_mb[snp_index_5B] > 110 &
#     wheat_data_phys$snp$pos_mb[snp_index_5B] < 210)
# ]
# snp_index_6A <- which(wheat_data_phys$snp$chrom == "6A")
# haplo_id_6A <- wheat_data_phys$snp$id[snp_index_6A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_6A] > 170 &
#     wheat_data_phys$snp$pos_mb[snp_index_6A] < 445)
# ]
# snp_index_6B <- which(wheat_data_phys$snp$chrom == "6B")
# haplo_id_6B <- wheat_data_phys$snp$id[snp_index_6B][
#   which(wheat_data_phys$snp$pos_mb[snp_index_6B] > 250 &
#     wheat_data_phys$snp$pos_mb[snp_index_6B] < 380)
# ]
# snp_index_7A <- which(wheat_data_phys$snp$chrom == "7A")
# haplo_id_7A <- wheat_data_phys$snp$id[snp_index_7A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_7A] > 310 &
#     wheat_data_phys$snp$pos_mb[snp_index_7A] < 445)
# ]
# haplo_ids <- c(
#   haplo_id_1A, haplo_id_2A, haplo_id_4A, haplo_id_5B, haplo_id_6A, haplo_id_6B,
#   haplo_id_7A
# )
# haplo_index_snps_phys <- match(haplo_ids, wheat_data_phys$snp$id)
# haplo_index_snps_gen <- match(haplo_ids, wheat_data_gen$snp$id)

# # create a column in the snp_data set that contains D values for only
# # those markers in the extended haplotypes, NA for all else
# phys_gen_snp_pos <- phys_gen_snp_pos %>% mutate(haplo_phys = phys)
# phys_gen_snp_pos$haplo_phys[-haplo_index_snps_phys] <- NA
# phys_gen_snp_pos <- phys_gen_snp_pos %>% mutate(haplo_gen = gen)
# phys_gen_snp_pos$haplo_gen[-haplo_index_snps_gen] <- NA

      # geom_point(
      #   aes(haplo_phys, haplo_gen),
      #   colour = "black",
      #   shape = 1, size = 0.75
      # )