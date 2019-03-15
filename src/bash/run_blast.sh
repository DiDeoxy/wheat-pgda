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
# # download and unzip refseq and gtf
# # echo "download_unzip_refseq_and_gtf"
# mkdir refseq
# wget -P refseq/ ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
# gunzip refseq Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
# wget -P refseq/ ftp://ftp.ensemblgenomes.org/pub/plants/release-42/gtf/triticum_aestivum/Triticum_aestivum.IWGSC.42.gtf.gz
# gunzip refseq/Triticum_aestivum.IWGSC.42.gtf.gz


# # 2
# # extract genes from refseq delete reference
# echo "extract_genes_from_reference_and_delete"
# mkdir -p data/raw/blast_db
# python src/python/extract_gene_sequences.py \
#     refseq/Triticum_aestivum.IWGSC.42.gtf \
#     refseq/Triticum_aestivum.IWGSC.dna.toplevel.fa \
#     data/raw/blast_db/ref_seq_v1_genes.fasta
# rm -r refseq

# # 3
# # make blast database from cds
# echo "make_blast_db"
# makeblastdb -in data/raw/blast_db/ref_seq_v1_genes.fasta \
#     -dbtype nucl \
#     -parse_seqids

# 4
# concat all resi and pheno fasta
echo "concat_resi"
rm data/raw/gene_seqs/resi/all.fasta
cat data/raw/gene_seqs/resi/*.fasta > data/raw/gene_seqs/resi/all.fasta
echo "concat_pheno"
rm data/raw/gene_seqs/pheno/all.fasta
cat data/raw/gene_seqs/pheno/*.fasta > data/raw/gene_seqs/pheno/all.fasta

# 5
# align resi genes
echo "align_resi_genes"
mkdir -p data/intermediate/blast
blastn \
    -num_threads 4 \
    -query data/raw/gene_seqs/resi/all.fasta \
    -db data/raw/blast_db/ref_seq_v1_genes.fasta \
    -outfmt "6 qseqid sseqid bitscore pident evalue qlen length" \
    > data/intermediate/blast/resi.txt

# 6
# get top resi alignments
echo "get_top_resi_alignments"
python src/python/get_top_alignments.py \
    data/intermediate/blast/resi.txt \
    data/raw/blast_db/ref_seq_v1_genes.fasta \
    data/intermediate/blast/top_resi.csv

# 7
# align pheno genes
echo "align_pheno_genes"
blastn \
    -num_threads 4 \
    -query data/raw/gene_seqs/pheno/all.fasta \
    -db data/raw/blast_db/ref_seq_v1_genes.fasta \
    -outfmt "6 qseqid sseqid bitscore pident evalue qlen length" \
    > data/intermediate/blast/pheno.txt

# 8
# get top pheno alignments
echo "get_top_pheno_alignments"
python src/python/get_top_alignments.py \
    data/intermediate/blast/pheno.txt \
    data/raw/blast_db/ref_seq_v1_genes.fasta \
    data/intermediate/blast/top_pheno.csv

# 9
# convert genbank ids of aligned gene sequences to gene names
echo "convert_genbank_ids"
bash src/bash/convert_genbank_ids_to_gene_names.sh

###############################################################################
# gotta do things different since fhb1 is missing from CS

# # 1
# # download 3B chromsome
# echo "download_chromosome_3B"
# mkdir refseq
# # wget -P refseq ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.chromosome.3B.fa.gz
# # gunzip refseq/Triticum_aestivum.IWGSC.dna.chromosome.3B.fa.gz
# echo 'take first two million lines'
# head -n 2000000 refseq/Triticum_aestivum.IWGSC.dna.chromosome.3B.fa > \
#     refseq/Triticum_aestivum.IWGSC.dna.chromosome.3B_first_two_million.fa

# # 2
# # make blast database from cds
# echo "make_3B_blast_db"
# makeblastdb -in refseq/Triticum_aestivum.IWGSC.dna.chromosome.3B_first_two_million.fa \
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
#     -db refseq/Triticum_aestivum.IWGSC.dna.chromosome.3B_first_two_million.fa \
#     -outfmt "6 qseqid sseqid bitscore pident evalue qlen length sstart send qstart qend" \
#     > data/intermediate/blast/FHB1.txt
# # rm -r refseq

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