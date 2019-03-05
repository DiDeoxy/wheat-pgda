"""
Max H, Oct 15, 2017
This script formats the blast output
"""

from sys import argv
import pandas as pd


def read(blast_alignment):
    """Read the blast data and sort."""
    blast = pd.read_csv(blast_alignment, sep='\t', names=[
        'name', 'chr', 'bitscore', 'pident', 'evalue', 'sstart', 'send'
    ])
    blast['name'] = blast['name'].str[:-2]
    blast['bitscore'] = blast['bitscore'].astype('int32')
    blast['pos'] = ((blast['sstart'] + blast['send'])/2).astype('int32')
    blast = blast.drop(['sstart', 'send'], 1)
    return blast.sort_values('bitscore', ascending=False)


def main():
    """Initialize the logic of the program."""
    blast_file, out_file = argv[1], argv[2]

    with open(out_file, 'w') as best:
        read(blast_file).to_csv(best, index=False)


if __name__ == "__main__":
    main()
