################################################################################
# create gene sequences form gff and refse v1, create blast db, align genes,
# extract top alignments, convert ncbi identifiers to gene names

# # 0
# # install needed software
# echo "install_blast"
# wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.8.1+-2.x86_64.rpm
# sudo dnf install ncbi-blast-2.8.1+-2.x86_64.rpm
# rm ncbi-blast-2.8.1+-2.x86_64.rpm
# echo "setup_python_venv"
#
# python3 -m venv src/python/ALIGN
# source src/python/ALIGN/bin/activate
# python -m pip install --upgrade pip
# python -m pip install biopython pandas pyfaidx
# python -m pip install git+https://github.com/daler/gffutils.git
# deactivate

# # 1
# # build gtf db and extract gene sequences
# echo "extract_gene_sequences"
# mkdir -p data/raw/blast_db
# python src/python/extract_gene_sequences.py \
#     refseq/Triticum_aestivum.IWGSC.42.gtf \
#     refseq/Triticum_aestivum.IWGSC.dna.toplevel.fa \
#     data/raw/blast_db/ref_seq_v1_genes.fasta

# # 2
# # make blast database from cds
# echo "make_blast_db"
# makeblastdb -in data/raw/blast_db/ref_seq_v1_genes.fasta \
#     -dbtype nucl \
#     -parse_seqids

# # 3
# # align resi genes
# echo "align_resi_genes"
# mkdir -p data/intermediate/blast
# blastn \
#     -num_threads 4 \
#     -query data/raw/gene_seqs/resi/all.fasta \
#     -db data/raw/blast_db/ref_seq_v1_genes.fasta \
#     -outfmt "6 qseqid sseqid bitscore pident evalue" \
#     > data/intermediate/blast/resi.txt

# # 4
# # align rpheno genes
# echo "align_pheno_genes"
# blastn \
#     -num_threads 4 \
#     -query data/raw/gene_seqs/pheno/all.fasta \
#     -db data/raw/blast_db/ref_seq_v1_genes.fasta \
#     -outfmt "6 qseqid sseqid bitscore pident evalue" \
#     > data/intermediate/blast/pheno.txt

# # 5
# # get top pheno alignments
# echo "get_top_resi_alignments"
# python src/python/get_top_alignments.py \
#     results/blast/resi.txt \
#     data/raw/blast_db/ref_seq_v1_genes.fasta \
#     data/intermediate/blast/top_resi.csv

# # 6
# # get top resi alignments
# echo "get_top_pheno_alignments"
# python src/python/get_top_alignments.py \
#     results/blast/pheno.txt \
#     data/raw/blast_db/ref_seq_v1_genes.fasta \
#     data/intermediate/blast/top_pheno.csv

# # 7
# # convert genbank ids of aligned gene sequences to gene names
# echo "convert_genbank_ids"
# bash src/bash/convert_genbank_ids_to_gene_names.sh

# ################################################################################
# # format and process data R

# # 0
# # install needed R packages
# echo "install_R_packages"
# Rscript src/R/install_packages.R

# # 1
# # take the two genetic maps we have and compare them to identify markers
# # assigned to the same linkage groups
# echo "create_filtered_map"
# Rscript src/R/create_filtered_map.R

# # 2
# # create a tibble using the consensus markers from step 1 with the physical
# # positions from the gmap alignment of the probes supplied by the IWGSC
# # with the genotype information for those markers supplied by Dr. Pozniak
# echo "create_phys_map_genotypes_tibble"
# Rscript src/R/create_phys_map_genotypes_tibble.R

# # 3
# # create gds format files from the table created in step 2 one with gmap
# # postions, the other with pozniak genetic map positions
# echo "create_gds"
# Rscript src/R/create_gds.R

# # 4
# # create a series of gds files, first pruning sets of samples with >99% IBS
# # then removing markers with a MAF < 0.05 or a MR > 0.10, finally LD prune
# echo "sample_maf_mr_and_ld_prune_data"
# Rscript src/R/sample_maf_mr_and_ld_prune_data.R

# # 5
# # calculate HDBSCAN clusters
# echo "hdbscan_cluster"
# Rscript src/R/hdbscan_cluster.R

# # 6
# # create geninds for calculating Jost's Ds and AMOVA phis between groups
# echo "create_geninds"
# Rscript src/R/create_geninds.R

# # 7
# # calculate Jost's D values
# echo "calc_josts_d"
# Rscript src/R/calc_josts_d.R

# ################################################################################
# # create and output results

# # 1
# # calculate and plot statistics of the maf and mr pruned phys map and the ld
# # pruned map
# echo "calc_and_plot_map_stats"
# Rscript src/R/calc_and_plot_map_stats.R

# # 2
# # plot a UPGMA dendrogram with hdbscan and metadata in rows of rectangels
# # arranged around it
# echo "plot_hdbscan_upgma_dend"
# Rscript src/R/plot_hdbscan_upgma_dend.R

# # 3
# # calculate AMOVA phi values for various hierarchical AMOVAs
# echo "calc_sample_groupings_AMOVA_phi"
# Rscript src/R/calc_sample_groupings_AMOVA_phi.R

#4
# plot the josts D values of markers between the varieties of the major
# phenptypes of the three largest clusters
echo "plot_clustered_phenos_markers_josts_ds_with_genes"
Rscript src/R/plot_clustered_phenos_markers_josts_ds_with_genes.R

# # 5
# # plot marker physical postions against genetic postion coloured by eh value
# echo "plot_phys_vs_gen_pos_with_eh"
# Rscript src/R/plot_phys_vs_gen_pos_with_eh.R

# # 6
# # plot ld heatmaps of each chromosome
# echo "plot_LD_heatmaps"
# Rscript src/R/plot_LD_heatmaps.R