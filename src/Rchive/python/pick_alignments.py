"""
Max H, Oct 15, 2017
This script formats the blast output
"""

from sys import argv
import pandas as pd


def read(blast_alignment):
    """Read the blast data and sort."""
    blast = pd.read_csv(blast_alignment, sep='\t', names=[
        'name', 'chr', 'query_length', 'align_length', 'pident', 'evalue',
        'bitscore', 'sstart', 'send'])
    blast = blast.set_index(['name', 'bitscore']).sort_index(ascending=False)
    return blast


def top_blast(alignments):
    """Find the best alignments."""
    best_position = list()
    for name, group in alignments.groupby(level=0):
        if isinstance(group, pd.DataFrame):
            count = 0
            for row in group.itertuples():
                best_position.append(
                    ",".join(
                        [name[:-2], str(row[1]), str(int((row[6] + row[7])/2)),
                         str(row[2]), str(row[3]), str(row[4])]
                    )
                )
                if count == 4:
                    break
                count = count + 1
    return best_position


def main():
    """Initialize the logic of the program."""
    blast_file, out_file = argv[1], argv[2]

    best_position = top_blast(read(blast_file))

    with open(out_file, 'w') as best:
        for position in best_position:
            best.write("%s\n" % position)


if __name__ == "__main__":
    main()
