#!/usr/bin/env python

#######################################################################
#
#   contigFilter.py
#   2015 James Stapleton
#
#   This script parses a list of synthetic reads,
#       and writes those those longer than minlength and
#       shorter than maxlength to a new file
#
#######################################################################

import argparse


def main(contig_file, minlength, maxlength, outfile):

    with open(contig_file, 'r') as contiglist, open (outfile, 'w') as filtered:
        trigger = 0
        contig = []
        for line in contiglist:
            #if line is a barcode, save it
            if line[:4] == '>Bar':
                barcode = line
            # If line is blank, end of barcode.
            # If triggered, write contig and start over.
            elif line[0] == '\n':
                if trigger:
                    if barcode:
                        filtered.write(barcode)
                    for stuff in contig:
                        filtered.write(stuff)
                    filtered.write('\n')
                trigger = 0
                contig = []
            # If line is a contig head, check length.
            # If big, trigger and start contig
            elif line[:5] == '>NODE':
                if trigger:
                    if barcode:
                        filtered.write(barcode)
                        barcode = 0
                    for stuff in contig:
                        filtered.write(stuff)
                    trigger = 0
                contig = []
                bases = int(line.split('_')[3])
                if bases > int(minlength) and bases < int(maxlength):
                    trigger = 1
                    contig.append(line)
            #if line is sequence, append to growing contig
            else:
                if trigger:
                    contig.append(line)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('contig_file')
    parser.add_argument('--minlength', default=0,
            help='minimum contig length to be written')
    parser.add_argument('--maxlength', default=99999999999999999999,
            help='maximum contig length to be written')
    parser.add_argument('--outfile', default='filtered_contigs.txt',
            help='output file name')
    args = parser.parse_args()

    main(args.contig_file, args.minlength, args.maxlength, args.outfile)
