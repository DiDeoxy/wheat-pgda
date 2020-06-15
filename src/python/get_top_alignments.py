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
    alignments = pd.read_csv(
        blast_file, sep='\t',
        names=['query', 'gene', 'bitscore', 'pident', 'evalue',
               'query_length', 'align_length']
    )
    alignments['query'] = alignments['query'].str[:-2]
    alignments['bitscore'] = alignments['bitscore'].astype('int32')
    return alignments.set_index(
        ['query', 'bitscore']
    ).sort_index(ascending=False).reset_index('bitscore')


def top_blast(alignments):
    """Find the best alignments."""
    top_aligns = list()
    for query, group in alignments.groupby('query'):
        top_aligns.append(
            group.loc[
                (group['pident'] >= 90) & (group['bitscore'] >= 800) &
                (group['bitscore'] >= group['bitscore'].iloc[0] * 0.5)
            ]
        )
    top_aligns = pd.concat(top_aligns).reset_index().set_index(
        ['gene', 'bitscore']).sort_index(ascending=False)
    top_aligns_by_gene = list()
    for gene, group in top_aligns.groupby('gene'):
        top_aligns_by_gene.append(group.iloc[0:1])
    return pd.concat(top_aligns_by_gene).reset_index()


def combine(alignments, genes):
    """Parse the fasta file"""
    loci = list()
    for gene in alignments['gene']:
        locus = genes[gene].description.split()
        loci.append([locus[1], int((int(locus[2]) + int(locus[3])) / 2)])
    return alignments.reset_index().merge(
        pd.DataFrame(loci, columns=['chr', 'pos']), left_index=True,
        right_index=True
    )[['query', 'gene', 'chr', 'pos', 'bitscore', 'pident', 'evalue',
       'query_length', 'align_length']]


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
