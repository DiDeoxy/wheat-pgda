"""
Max H, Jun 8, 2017
This script formats the blast output
"""

from sys import argv
import pandas as pd

def parse_blast(file):
    """this function parses the blast output"""
    with open(file) as alignments:
        blast = list()
        for alignment in alignments:
            (name, chromo_part, query_length, align_length,
             pident, evalue, bitscore, sstart, send) = alignment.split()
            if chromo_part == "chrUn":
                continue
            else:
                (chromo, part) = chromo_part.split("_")
                chromo = chromo[-2:]
            name = name[:-2]
            blast.append(",".join([name, chromo, part, query_length, align_length,
                                   pident, evalue, bitscore, sstart, send]))
    return pd.DataFrame(blast, index_col=0, names=[
        'name', 'chr', 'part', 'query_length',
        'align_length', 'pident', 'evalue',
        'bitscore', 'sstart', 'send'])

def seek(blast):
    """this function compares the blast output with the gen map"""

    print(blast.index)

def main():
    """This runs the program"""
    blast_file, out_file = argv[1], argv[2]

    blast = parse_blast(blast_file)

    seek(blast)

    with open(out_file, 'w') as formatted:
        for alignment in blast:
            formatted.write("%s\n" % alignment)

if __name__ == "__main__":
    main()
