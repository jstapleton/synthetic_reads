#!/usr/bin/env python

################################################################
#
#   unique_barcodes.py
#   2015 James Stapleton
#
#   Compare all barcodes in a contig list against all others.
#   If two barcodes are close, discard the shorter one
#   Before running this, use longest_contig_writer.py
#   to make sure there is only 1 contig/barcode
#
################################################################

import argparse


def hamdist(str1, str2):
    """Count the # of differences between equal-length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs


def main(contig_file, outfile):

# 1. Make a list of lists containing barcode and synthetic read length
    barcode_length_pairs = []
    with open(contig_file, 'r') as contigs:
        contiglist = contigs.readlines()
        for line in contiglist:
            #if line is a barcode, save it
            if line[:4] == ">Bar":
                barcode = line[10:26]
            #if line is a contig head, append barcode and length
            elif line[:5] == '>NODE':
                bases = int(line.split('_')[3])
                barcode_and_length = [barcode, bases]
                barcode_length_pairs.append(barcode_and_length)

# 2. Compare each barcode against the others
#    When two are close, write the one with the shorter read to a list
    bad_barcodes = set()
    for i in range(len(barcode_length_pairs)-1):
        for j in range(i+1, len(barcode_length_pairs)):
            hamming = hamdist(barcode_length_pairs[i][0], barcode_length_pairs[j][0])
            if hamming < 3 and hamming != 0:
                if barcode_length_pairs[i][1] > barcode_length_pairs[j][1]:
                    bad_barcodes.add(barcode_length_pairs[j][0])
                else:
                    bad_barcodes.add(barcode_length_pairs[i][0])
                    break
    print "found", str(len(bad_barcodes)), "bad barcodes out of", \
        str(len(barcode_length_pairs)), "total barcodes"

# 3. Remove the duplicate barcodes from the contig_list.
# Go through contiglist looking for barcodes on the list,
# printing everything else to outfile
    trigger = 1  # if 0, we are in a barcode to be removed
    with open(contig_file, 'r') as contigs:
        with open(outfile, 'w') as outf:
            for line in contigs:
                if line[:2] == '>B':
                    barcode = line[10:26]
                    if barcode in bad_barcodes:
                        trigger = 0
                    else:
                        trigger = 1
                if trigger:
                    outf.write(line)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('contig_file')
    parser.add_argument('--outfile', default='unique_contig_list.txt')
    args = parser.parse_args()

    main(args.contig_file, args.outfile)
