#!/usr/bin/env python

################################################################
#
#   longest_contig_writer.py
#   2015 James Stapleton
#
#   Takes a file of synthetic reads.
#   Writes a new file with the single longest synthetic read
#   under each barcode group.
#
################################################################

import argparse


def main(contig_file, outfile):
    with open(contig_file, 'r') as contigs, open(outfile, 'a') as maxlist:
        max = 0
        max_contig = []
        trigger = 0
        for line in contigs:
            #if line is a barcode, save it
            if line[:4] == ">Bar":
                barcode = line
            #if line is a contig head, compare length to max
            elif line[:5] == '>NODE':
                bases = int(line.split('_')[3])
                if bases > max:
                    trigger = 1
                    max = bases
                    max_contig = []
                    max_contig.append(line)
                else:
                    trigger = 0
            # if line is blank, end of barcode:
            # write out barcode and longest contig
            elif line[0] == '\n':
                maxlist.write(barcode)
                for line in max_contig:
                    maxlist.write(line)
                maxlist.write('\n')
                # reset
                max = 0
                max_contig = []
                trigger = 0
            else:
                if trigger:
                    max_contig.append(line)
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('contig_file')
    parser.add_argument('--outfile', default='longest_contigs.txt')
    args = parser.parse_args()

    main(args.contig_file, args.outfile)
