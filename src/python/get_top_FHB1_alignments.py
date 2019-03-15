"""
Max H, Oct 15, 2017
This script formats the blast output
"""

from sys import argv
import pandas as pd


def parse_blast(blast_file):
    """Read the blast data and sort."""
    alignments = pd.read_csv(
        blast_file, sep='\t',
        names=['query', 'chr', 'bitscore', 'pident', 'evalue',
               'query_length', 'align_length', 'chr_start', 'chr_end',
               'query_start', 'query_end']
    )
    alignments['query'] = alignments['query'].str[:-2]
    alignments['bitscore'] = alignments['bitscore'].astype('int32')
    alignments['gene'] = 'None'
    alignments['pos'] = (
        (alignments['chr_start'] + alignments['chr_end'])/2).astype('int32')
    return alignments.set_index(
        ['query', 'bitscore']
    ).sort_index(ascending=False).reset_index('bitscore')


def top_alignments(alignments):
    """Find the best alignments."""
    top_aligns = list()
    for query, group in alignments.groupby('query'):
        if isinstance(group, pd.DataFrame):
            top_aligns.append(
                group.loc[group['bitscore'] >= group['bitscore'].iloc[0] * 0.5]
            )
    return pd.concat(top_aligns).reset_index()[
        ['query', 'gene', 'chr', 'pos', 'bitscore', 'pident', 'evalue',
         'query_length', 'align_length', 'chr_start', 'chr_end', 'query_start',
         'query_end']
    ]


def main():
    """Initialize the logic of the program."""
    blast_file, out_file = argv[1], argv[2]

    with open(out_file, 'w') as best:
        top_alignments(parse_blast(blast_file)).to_csv(best, index=False)


if __name__ == "__main__":
    main()
