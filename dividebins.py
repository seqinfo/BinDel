import argparse
import sys

# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com

parser = argparse.ArgumentParser(description='Generate bins.')

parser.add_argument('--infile',
                    type=argparse.FileType('r', encoding="UTF-8"),
                    help="File in the TSV format (UTF-8) with columns: 'chr', 'start', 'end', 'focus', 'length'.")

parser.add_argument('--outfile',
                    nargs="?",
                    type=argparse.FileType('w', encoding="UTF-8"),
                    default=sys.stdout,
                    help="Binned output file.")

args = parser.parse_args()

# Skip header
args.infile.readline()
header = "\t".join(['chr', 'start', 'end', 'focus'])

args.outfile.write(header.strip() + "\n")
for line in args.infile:
    line = line.strip().split("\t")
    chromosome, start, end, focus, bin_width = line[0], int(line[1]), int(line[2]), line[3], int(line[4])

    while start + bin_width < end:
        args.outfile.write(f"{chromosome}\t{start}\t{(min(start + bin_width - 1, end))}\t{focus}\n")
        start = start + bin_width

args.infile.close()
args.outfile.close()
