# make the blast db
makeblastdb -in data/raw/blastdb/Triticum_aestivum.IWGSC.cds.all.fa \
    -dbtype nucl \
    -parse_seqids