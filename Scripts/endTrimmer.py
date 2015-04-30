#!/usr/bin/env python

############################################################################
#
#   endTrimmer.py
#   2015 James Stapleton
#
#   Trims a specified number of nucleotides from each end of each
#       synthetic long read in a list, prints a new file
#
############################################################################

import argparse


def main(contig_file, trim_length, outfile):
    TRIM_LENGTH = int(trim_length)
    with open(contig_file, 'r') as contigs:
        with open(outfile, 'w') as outfile:
            sequence = ''
            for line in contigs:
                #if line is a barcode, print it
                if line[:4] == ">Bar":
                    outfile.write(line)
                #if line is a contig head, print it,
                # write out trimmed sequence
                elif line[:5] == '>NODE':
                    if sequence:
                        outfile.write(sequence[TRIM_LENGTH:-TRIM_LENGTH])
                        outfile.write('\n')
                        sequence = ''
                    outfile.write(line)
                # if line is blank, end of barcode:
                # write out trimmed sequence
                elif line[0] == '\n':
                    outfile.write(sequence[TRIM_LENGTH:-TRIM_LENGTH])
                    outfile.write('\n\n')
                    sequence = ''
                else:
                    sequence = sequence + line.rstrip()

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('contig_file')
    parser.add_argument('trim_length')
    parser.add_argument('--outfile', default='trimmed_contigs.txt')
    args = parser.parse_args()

    main(args.contig_file, args.trim_length, args.outfile)
