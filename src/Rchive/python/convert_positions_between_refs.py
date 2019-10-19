"""Convert between chrPart and chr positions for the wheat reference genome."""

import argparse
import csv
import re


def main():
    """Initialize the program."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "ref_pos",
        help=('file containing the conversion between chrPart and chr '
              'positions'))
    parser.add_argument(
        "align_pos",
        help="file containing the gene alignments in csv format")
    parser.add_argument("outfile", help="the file to print the conversion to")
    args = parser.parse_args()

    ref_pos = {}
    with open(args.ref_pos) as csvfile:
        contents = csv.reader(csvfile, delimiter='\t')
        for row in contents:
            row[0] = re.sub("chr", "", row[0])
            row[4] = re.sub("chr", "", row[4])
            ref_pos[row[0]] = row[1:]

    align_pos = {}
    with open(args.align_pos) as csvfile:
        contents = csv.reader(csvfile, delimiter=',')
        next(contents)
        for row in contents:
            align_pos[row[0]] = row[1:]

    with open(args.outfile, 'w', newline='') as csvfile:
        csvfile.write("Name,ID,Chr,Position,seq_length,align_length,%id\n")
        contents = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_ALL)
        for k, v, in align_pos.items():
            v[2] = str(int(ref_pos[v[1]][4]) + int(v[2]))
            # if int(v[2]) > int(ref_pos[v[1]][5]):
            #     print("bad")
            v[1] = ref_pos[v[1]][3]
            v.insert(0, k)
            contents.writerow(v)


if __name__ == "__main__":
    main()
