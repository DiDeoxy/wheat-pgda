source(file.path("repos", "wheat-pgda", "src", "R", "file_paths.R"))
library("biomaRt")
import::from(dplyr, "arrange", "mutate", "select")
import::from(GenomicFeatures, "genes", "makeTxDbFromGFF")
import::from(GenomicRanges, "findOverlaps", "GRanges", "mcols", "mcols<-")
import::from(IRanges, "extractList", "IRanges")
import::from(magrittr, "%>%", "%<>%")
import::from(pgda, "max_lengths", "snpgds_parse", "load_genes")
import::from(plyr, "rbind.fill")
import::from(readr, "read_rds")
import::from(refGenome, "getGenePositions", "ensemblGenome", "read.gtf")
import::from(S4Vectors, "queryHits", "unstrsplit")
library("topGO")
import::from(tibble, "as_tibble")

# load the data from the gds object
phys_data <- snpgds_parse(phys_gds)

# load some data
josts_d <- read_rds(josts_d)

# get chomr names order
chroms <- as.vector(
  t(outer(as.character(1:7), c("A", "B", "D"), paste, sep = ""))
)

# construct the ranges to extract genes from
sig_markers <- phys_data$snp[which(josts_d > quantile(josts_d, 0.80)), ]
sig_ranges <- GRanges(
    sig_markers$chrom,
    IRanges(sig_markers$pos - 500000, sig_markers$pos + 500000)
)

# gene_ranges <- makeTxDbFromGFF("/workspace/repos/wheat-pgda/src/bash/ref_seq/Triticum_aestivum.IWGSC.42.gtf") %>% GenomicFeatures::genes()

overlaps <- findOverlaps(sig_ranges, gene_ranges, ignore.strand = TRUE)

all_genes <- mcols(gene_ranges)$gene_id

sig_genes <- extractList(all_genes, as(overlaps, "List")) %>% unlist()

gene_list <- factor(as.integer(all_genes %in% sig_genes))
names(gene_list) <- all_genes

gene_go_id <- getBM(
  attributes = c("ensembl_gene_id", "go_id"),
  mart = useMart(
    "plants_mart",
    dataset = "taestivum_eg_gene",
    host = "plants.ensembl.org"
  )
) %>%
  .[.$go_id != '', ]

gene_go_id_list <- by(gene_go_id$go_id, gene_go_id$ensembl_gene_id, as.character)

go_data <- new(
  "topGOdata", ontology = "BP", allGenes = gene_list,
  annot = annFUN.gene2GO, gene2GO = gene_go_id_list
)

result_fisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")

all_res <- GenTable(
  go_data, classicFisher = result_fisher,
  orderBy = "result_fisher", ranksOf = "classicFisher",
  topNodes = length(usedGO(object = go_data))
)
all_res$classicFisher <- p.adjust(all_res$classicFisher, method = "bonferroni")
head(all_res, n = 10)
printGraph(
  go_data, result_fisher, firstSigNodes = 6, fn.prefix = "fisher",
  useInfo = "all", pdfSW = TRUE
)

genes_auxin <- genesInTerm(go_data, "GO:0009733") 
sig_genes_auxin <- genes_auxin[[1]][genes_auxin[[1]] %in% sigGenes(go_data)]

genes_wounding <- genesInTerm(go_data, "GO:0009611")
sig_genes_wounding <- genes_wounding[[1]][genes_wounding[[1]] %in% sigGenes(go_data)]

# genes_neg_reg_endo <- genesInTerm(go_data, "GO:0010951")
# sig_genes_neg_reg_endo <- genes_neg_reg_endo[[1]][genes_neg_reg_endo[[1]] %in% sigGenes(go_data)]

# genes_biosynth <- genesInTerm(go_data, "GO:0070681")
# sig_genes_biosynth <- genes_biosynth[[1]][genes_biosynth[[1]] %in% sigGenes(go_data)]

sig_genes <- rbind(
  gene_ranges[which(mcols(gene_ranges)$gene_id %in% sig_genes_auxin), ] %>%
    as_tibble() %>%
    mutate(type = "Response to Auxin", base = 0.5, pos_mb = (start + end) / 2e6) %>%
    dplyr::select(id = gene_id, chrom = seqnames, base, type, pos_mb),
  gene_ranges[which(mcols(gene_ranges)$gene_id %in% sig_genes_wounding), ] %>%
    as_tibble() %>%
    mutate(type = "Response to Wounding", base = 0.5, pos_mb = (start + end) / 2e6) %>%
    dplyr::select(id = gene_id, chrom = seqnames, base, type, pos_mb)
) %>% arrange(chrom, pos_mb)

unique(sig_genes$chrom)

write.table(
  sig_genes,
  file.path("/workspace", "data", "intermediate", "go_genes.csv"), sep = "," , row.names = FALSE
)

# wounding_genes <- gene_ranges[which(mcols(gene_ranges)$gene_id %in% sig_genes_wounding), ] %>%
#   as_tibble() %>%
#   mutate(type = "Response to Wounding", base = 0.5, pos_mb = (start + end) / 2e6) %>%
#   dplyr::select(id = gene_id, chrom = seqnames, base, type, pos_mb)

# neg_reg_endo_genes <- gene_ranges[which(mcols(gene_ranges)$gene_id %in% sig_genes_neg_reg_endo), ] %>%
#   as_tibble() %>%
#   mutate(type = "Negative Rgulation Endopeptidase Activity", base = 0.5, pos_mb = (start + end) / 2e6) %>%
#   dplyr::select(id = gene_id, chrom = seqnames, base, type, pos_mb)

write.table(
  sig_genes %>%
    dplyr::select(`Gene Model` = id, Chromosome = chrom, `Pos in Mb` = pos_mb, Type = type) %>%
    arrange(Chromosome, `Pos in Mb`),
  file.path("results", "significant_genes.csv"), sep = "," , row.names = FALSE
)