import pandas as pd
from Bio import SeqIO
from sys import argv
"""
Max H, March 2, 2019
This script formats the blast output
modified from pick_alignments.py
"""


def parse_blast(blast_file):
    """Read the blast data and sort."""
    blast = pd.read_csv(
        blast_file, sep='\t',
        names=['query', 'subject', 'bitscore', 'pident', 'evalue']
    )
    blast['query'] = blast['query'].str[:-2]
    blast['bitscore'] = blast['bitscore'].astype('int32')
    blast = blast.set_index(
        ['query', 'subject', 'bitscore']
    ).sort_index(ascending=False)
    return blast


def combine(blast_data, fasta_data):
    """Parse the fasta file"""
    loci = list()
    for subject in blast_data['subject']:
        locus = fasta_data[subject].description.split()
        loci.append([locus[1], int((int(locus[2]) + int(locus[3])) / 2)])
    combined = blast_data.reset_index().merge(
        pd.DataFrame(loci, columns=['chr', 'pos']), left_index=True,
        right_index=True
    )[['query', 'subject', 'chr', 'pos', 'bitscore', 'pident', 'evalue']]
    combined['subject'] = combined['subject'].str[:-2]
    return combined


def top_blast(alignments):
    """Find the best alignments."""
    top_by_q_s = list()
    for q_s, group in alignments.groupby(level=['query', 'subject']):
        if isinstance(group, pd.DataFrame):
            top_by_q_s.append(group.iloc[0:1, ])
    top_by_q_s = pd.concat(top_by_q_s).reset_index(
        level='subject').sort_index(ascending=False)

    top_aligns = list()
    for query, group in top_by_q_s.groupby(level=['query']):
        if isinstance(group, pd.DataFrame):
            top_aligns.append(group.iloc[0:5, ])
    return pd.concat(top_aligns)


def main():
    """Initialize the logic of the program."""
    blast_file, fasta_file, out_file = argv[1], argv[2], argv[3]

    with open(out_file, 'w') as best:
        combine(
            top_blast(parse_blast(blast_file)),
            SeqIO.index(fasta_file, 'fasta')
        ).to_csv(best, index=False)


if __name__ == "__main__":
    main()
