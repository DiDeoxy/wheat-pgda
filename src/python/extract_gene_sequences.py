import gffutils
import pyfaidx
import os
from sys import argv
"""
Max H, March 4, 2019
This script extracts genes from the wheat refseq v1.0
"""


def main():
    """Initialize the logic of the program."""
    gtf_file, fasta_file, out_file = argv[1], argv[2], argv[3]

    gtf_db = gtf_file + "_gffutils.db"

    if not os.path.isfile(gtf_db):
        gffutils.create_db(
            gtf_file, dbfn=gtf_db, disable_infer_genes=True,
            disable_infer_transcripts=True
        )

    db = gffutils.FeatureDB(gtf_db)
    fasta = pyfaidx.Fasta(fasta_file)

    genes = list()
    for gene in db.features_of_type("gene"):
        genes.append(
            " ".join([">" + gene.id, gene.chrom, str(gene.start),
                      str(gene.end), "\n" + gene.sequence(fasta)]
                     )
        )
    with open(out_file, 'w') as out:
        for gene in genes:
            out.write("%s\n" % gene)


if __name__ == "__main__":
    main()
