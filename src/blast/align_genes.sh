blastn \
    -num_threads 4 \
    -evalue 1E-35 \
    -query data/raw/gene_seqs/resi/all.fasta \
    -db data/raw/blastdb/Triticum_aestivum.IWGSC.cds.all.fa \
    -outfmt "6 qseqid sseqid bitscore pident evalue" \
    > results/blast/resi.txt

blastn \
    -num_threads 4 \
    -evalue 1E-35 \
    -query data/raw/gene_seqs/pheno/all.fasta \
    -db data/raw/blastdb/Triticum_aestivum.IWGSC.cds.all.fa \
    -outfmt "6 qseqid sseqid bitscore pident evalue" \
    > results/blast/pheno.txt