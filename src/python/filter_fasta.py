# import pandas as pd
from sys import argv
from Bio import SeqIO
import re
"""
Max H, March 16, 2019
THis script extracts fasta records with a ceratin pattern in their descriptor
"""


def main():
    """Initialize the logic of the program."""
    in_fasta, out_fasta, pattern = argv[1], argv[2], argv[3]

    regex = re.compile(pattern, re.IGNORECASE)

    seq_records = dict()
    for seq_record in SeqIO.parse(in_fasta, 'fasta'):
        if seq_record.id not in seq_records and regex.search(
                seq_record.description) and len(seq_record.seq) < 25000:
            seq_records[seq_record.id] = seq_record

    SeqIO.write(
        (seq_records[seq_record] for seq_record in seq_records), out_fasta,
        'fasta'
    )


if __name__ == "__main__":
    main()
