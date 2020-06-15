# create gene sequences form gff and refse v1, create blast db, align genes,
# extract top alignments, convert ncbi identifiers to gene names

# # 0
# # install needed software
# echo "install_blast"
# sudo apt-get install ncbi-blast+

# echo "setup_python_venv"
# sudo apt-get install -y python3-venv
# python3 -m venv src/python/ALIGN
source src/python/ALIGN/bin/activate
# python -m pip install --upgrade pip
# python -m pip install biopython pandas pyfaidx
# python -m pip install git+https://github.com/daler/gffutils.git

# # 1
# # download and unzip ref_seq and gtf
# # echo "download_unzip_ref_seq_and_gtf"
# mkdir ref_seq
wget -P ref_seq/ ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
gunzip ref_seq/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
wget -P ref_seq/ ftp://ftp.ensemblgenomes.org/pub/plants/release-42/gtf/triticum_aestivum/Triticum_aestivum.IWGSC.42.gtf.gz
gunzip ref_seq/Triticum_aestivum.IWGSC.42.gtf.gz

# # 2
# extract genes from ref_seq delete reference
echo "extract_genes_from_reference_and_delete"
mkdir -p ref_seq/blast_db
python src/python/extract_gene_sequences.py \
    ref_seq/Triticum_aestivum.IWGSC.42.gtf \
    ref_seq/Triticum_aestivum.IWGSC.dna.toplevel.fa \
    ref_seq/ref_seq_v1_genes.fasta
rm -r ref_seq

python .//extract_gene_sequences.py \
    Triticum_aestivum.IWGSC.42.gtf \
    Triticum_aestivum.IWGSC.dna.toplevel.fa \
    ref_seq_v1_genes.fasta

# # 3a
# # make blast database from whole genome
# echo "make_blast_db"
# makeblastdb -in ref_seq/Triticum_aestivum.IWGSC.dna.toplevel.fa \
#     -dbtype nucl \
#     -parse_seqids

# # 3b
# # make blast database from cds
# echo "make_blast_db"
# makeblastdb -in ref_seq/ref_seq_v1_genes.fasta \
#     -dbtype nucl \
#     -parse_seqids

# # 4a
# # concat resi
# echo "concat_ugt"
# rm data/raw/gene_seqs/resi/UGT/UGT.fasta
# cat data/raw/gene_seqs/resi/UGT/*.fasta > data/raw/gene_seqs/resi/UGT/UGT.fasta
# echo "filter_ugt"
# python src/python/filter_fasta.py \
#     data/raw/gene_seqs/resi/UGT/UGT.fasta \
#     data/raw/gene_seqs/resi/UGT_filtered.fasta \
#     "ugt|glycolsy|udp"
# echo "concat_resi"
# rm data/raw/gene_seqs/resi/all.fasta
# cat data/raw/gene_seqs/resi/*.fasta > data/raw/gene_seqs/resi/all.fasta

# # 5a
# # align resi genes
# echo "align_resi_genes"
# mkdir -p data/intermediate/blast
# blastn \
#     -num_threads 4 \
#     -query data/raw/gene_seqs/resi/all.fasta \
#     -db ref_seq/ref_seq_v1_genes.fasta \
#     -outfmt "6 qseqid sseqid bitscore pident evalue qlen length" \
#     > data/intermediate/blast/resi.txt

# # 6a
# # get top resi alignments
# echo "get_top_resi_alignments"
# python src/python/get_top_alignments.py \
#     data/intermediate/blast/resi.txt \
#     ref_seq/ref_seq_v1_genes.fasta \
#     data/intermediate/blast/top_resi.csv

# # 4b
# # concate pheno
# echo "concat_GLU"
# rm data/raw/gene_seqs/pheno/glu/glu.fasta
# cat data/raw/gene_seqs/pheno/glu/*.fasta > data/raw/gene_seqs/pheno/glu/glu.fasta
# echo "filter_GLU"
# python src/python/filter_fasta.py \
#     data/raw/gene_seqs/pheno/glu/glu.fasta \
#     data/raw/gene_seqs/pheno/GLU_filtered.fasta \
#     "glu|mw|weight"
# echo "filter_Pins"
# python src/python/filter_fasta.py \
#     data/raw/gene_seqs/pheno/Pins/pins.fasta \
#     data/raw/gene_seqs/pheno/Pins_filtered.fasta \
#     "pin"
# echo "concat_pheno"
# rm data/raw/gene_seqs/pheno/all.fasta
# cat data/raw/gene_seqs/pheno/*.fasta > data/raw/gene_seqs/pheno/all.fasta

# # 5b
# # align pheno genes
# echo "align_pheno_genes"
# blastn \
#     -num_threads 4 \
#     -query data/raw/gene_seqs/pheno/all.fasta \
#     -db ref_seq/ref_seq_v1_genes.fasta \
#     -outfmt "6 qseqid sseqid bitscore pident evalue qlen length" \
#     > data/intermediate/blast/pheno.txt

# # 6b
# # get top pheno alignments
# echo "get_top_pheno_alignments"
# python src/python/get_top_alignments.py \
#     data/intermediate/blast/pheno.txt \
#     ref_seq/ref_seq_v1_genes.fasta \
#     data/intermediate/blast/top_pheno.csv

# # 7
# # convert genbank ids of aligned gene sequences to gene names
# echo "convert_genbank_ids"
# bash src/bash/convert_genbank_ids_to_gene_names.sh

###############################################################################
# gotta do things different since fhb1 is missing from CS

# # 1
# # download 3B chromsome
# echo "download_chromosome_3B"
# mkdir ref_seq
# # wget -P ref_seq ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.chromosome.3B.fa.gz
# # gunzip ref_seq/Triticum_aestivum.IWGSC.dna.chromosome.3B.fa.gz
# echo 'take first two million lines'
# head -n 2000000 ref_seq/Triticum_aestivum.IWGSC.dna.chromosome.3B.fa > \
#     ref_seq/Triticum_aestivum.IWGSC.dna.chromosome.3B_first_two_million.fa

# # 2
# # make blast database from cds
# echo "make_3B_blast_db"
# makeblastdb -in ref_seq/Triticum_aestivum.IWGSC.dna.chromosome.3B_first_two_million.fa \
#     -dbtype nucl \
#     -parse_seqids

# # 3
# # align FHB1
# echo "align_FHB1"
# mkdir -p data/intermediate/blast
# blastn \
#     -num_threads 4 \
#     -max_hsps 10 \
#     -query data/raw/gene_seqs/resi/FHB1/all.fasta \
#     -db ref_seq/Triticum_aestivum.IWGSC.dna.chromosome.3B_first_two_million.fa \
#     -outfmt "6 qseqid sseqid bitscore pident evalue qlen length sstart send qstart qend" \
#     > data/intermediate/blast/FHB1.txt
# # rm -r ref_seq

# # 4
# # get top FHB1 alignment
# echo "get_top_FHB1_alignment"
# python src/python/get_top_FHB1_alignments.py \
#     data/intermediate/blast/FHB1.txt \
#     data/intermediate/blast/top_FHB1.csv

# # 5
# # convert the genbank ids to gene names
# echo "convert_genbank_ids"
# bash src/bash/convert_genbank_ids_to_gene_names.sh

###############################################################################
# 4AL est alignments

# # 5c
# # align 4AL13 genes
# echo "align_4AL13_genes"
# blastn \
#     -num_threads 4 \
#     -query data/raw/ests/gene_seqs/4AL13.fasta \
#     -db ref_seq/ref_seq_v1_genes.fasta \
#     -outfmt "6 qseqid sseqid bitscore pident evalue qlen length" \
#     > data/intermediate/blast/4AL13.txt

# # 6c
# # get top 4AL13 alignments
# echo "get_top_4AL13_alignments"
# python src/python/get_top_alignments.py \
#     data/intermediate/blast/4AL13.txt \
#     ref_seq/ref_seq_v1_genes.fasta \
#     data/intermediate/blast/top_4AL13.csv

# # 5d
# # align 4AL5 genes
# echo "align_4AL5_genes"
# blastn \
#     -num_threads 4 \
#     -query data/raw/ests/gene_seqs/4AL5.fasta \
#     -db ref_seq/ref_seq_v1_genes.fasta \
#     -outfmt "6 qseqid sseqid bitscore pident evalue qlen length" \
#     > data/intermediate/blast/4AL5.txt

# # 6d
# # get top 4AL5 alignments
# echo "get_top_4AL5_alignments"
# python src/python/get_top_alignments.py \
#     data/intermediate/blast/4AL5.txt \
#     ref_seq/ref_seq_v1_genes.fasta \
#     data/intermediate/blast/top_4AL5.csv

# # 5e
# # align 4AL4 genes
# echo "align_4AL4_genes"
# blastn \
#     -num_threads 4 \
#     -query data/raw/ests/gene_seqs/4AL4.fasta \
#     -db ref_seq/ref_seq_v1_genes.fasta \
#     -outfmt "6 qseqid sseqid bitscore pident evalue qlen length" \
#     > data/intermediate/blast/4AL4.txt

# # 6e
# # get top 4AL4 alignments
# echo "get_top_4AL4_alignments"
# python src/python/get_top_alignments.py \
#     data/intermediate/blast/4AL4.txt \
#     ref_seq/ref_seq_v1_genes.fasta \
#     data/intermediate/blast/top_4AL4.csv

##############################################################################
# Cereba alignments

# 4a
# concat resi
echo "concat_cereba"
rm data/raw/gene_seqs/cereba/cerebas.fasta
cat data/raw/gene_seqs/cereba/*.fasta > data/raw/gene_seqs/cereba/cerebas.fasta

# 5c
# align cereba genes
echo "align_cereba_genes"
blastn \
    -num_threads 4 \
    -query data/raw/gene_seqs/cereba/cerebas.fasta \
    -db ref_seq/Triticum_aestivum.IWGSC.dna.toplevel.fa \
    -outfmt "6 qseqid sseqid bitscore pident evalue qlen length sstart send" \
    > cereba_blast/cereba.txt

# # 6c
# # get top 4AL13 alignments
# echo "get_top_4AL13_alignments"
# python src/python/get_top_alignments.py \
#     data/intermediate/blast/4AL13.txt \
#     ref_seq/ref_seq_v1_genes.fasta \
#     data/intermediate/blast/top_4AL13.csv