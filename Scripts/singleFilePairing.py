#!/usr/bin/env python

###################################################################
#
#   singleFilePairing.py
#   2015 James Stapleton
#
#   For single-file barcode pairing:
#       Find paired barcodes from unfragmented molecules
#       within the same sample,
#       scan raw forward fastq for the correct pattern,
#       if found, print that sequence to a new file.
#
#   This is done separately from barcodeHasher.py because
#       this needs to be done before trimmomatic.
#
###################################################################

import argparse
import time
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def main(infile, afterBarcode, SCA2):
    BARCODE_LENGTH = 16
    COMPLEMENT_DICT = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N'}
    AFTER_BARCODE_SEQS = []
    AFTER_BARCODE_RC_SEQS = []
    if afterBarcode:
        for ABseq in afterBarcode.split(','):
            AFTER_BARCODE_SEQS.append(ABseq)
            AFTER_BARCODE_RC_SEQS.append(''.join([COMPLEMENT_DICT[base]
                                         for base in ABseq])[::-1])
            print "After barcode:", ABseq
    else:
        AFTER_BARCODE_SEQS.append('')
        AFTER_BARCODE_RC_SEQS.append('')
        print "No after-barcode sequence checking"

    if SCA2:
        LOSTUSEQS = ['AATTCCT', 'ATCGTTC']
    else:
        LOSTUSEQS = ['AGG', 'AATAGTT', 'ATGTGCATT']

    LOSTUSEQ_RCS = []
    for lostuseq in LOSTUSEQS:
        LOSTUSEQ_RCS.append(''.join([COMPLEMENT_DICT[base] for
                            base in lostuseq])[::-1])

    AFTER_BARCODES = 'AGATCGGAAGAGC'

    # initilize stuff
    readCount = 0
    numFound = 0
    start_time = time.time()

    with open(infile, 'rU') as fwdReads:
        # iterate over forward reads looking for paired barcodes
        with open("foundPairs.fastq", 'w') as outfile:
            for title, seq, qual in FastqGeneralIterator(fwdReads):
                # print time elapsed every 100000 reads processed
                readCount += 1
                if readCount % 100000 == 0:
                    print readCount
                    print time.time() - start_time, "seconds"
                    start_time = time.time()

                if len(seq) < BARCODE_LENGTH + BARCODE_LENGTH:
                    continue

                # pull out barcode1
                seq_F = seq[BARCODE_LENGTH:]

                # check for after_barcode_sequences
                if AFTER_BARCODE_SEQS:
                    test = after_barcode_seq(seq_F, AFTER_BARCODE_SEQS)
                    if test == 1:
                        continue
                    else:
                        seq_F = seq_F[len(test):]
                    if len(seq) < len(test) + BARCODE_LENGTH:
                        continue

                # check for lostU and lostUrc
                seq_F = trim_lost_U(seq_F, LOSTUSEQS)
                for lostuseq_RC in LOSTUSEQ_RCS[::-1]:
                    lostarray = [lostuseq_RC]
                    seq_F = trim_lost_U(seq_F, lostarray)
                # check for RC of after_barcode_sequences
                if AFTER_BARCODE_SEQS:
                    test = after_barcode_seq(seq_F, AFTER_BARCODE_RC_SEQS)
                    if test == 1:
                        continue
                    else:
                        seq_F = seq_F[len(test):]

                # pull out barcode2
                seq_F = seq_F[BARCODE_LENGTH:]

                # check after barcodes
                if seq_F[:len(AFTER_BARCODES)] == AFTER_BARCODES:
                    outfile.write('@'+title+'\n')
                    for record in seq, "+", qual:
                        outfile.write(record+'\n')
                    numFound += 1

    print str(readCount) + " total lines, " + str(numFound) + \
        " matching lines found"

    return 0


def after_barcode_seq(seq_F, AFTER_BARCODE_SEQS):
    "throw out sequences without the right sequence after the barcode"
    for ABseq in AFTER_BARCODE_SEQS:
        if seq_F[:len(ABseq)] == ABseq:
            return ABseq
    # if none of the allowed sequences are found,
    return 1


def trim_lost_U(seq_F, LOSTUSEQS):
    """ test for lost U at the 3' end of the PCR primer sequence """
    keepgoing = 1
    for lostuseq in LOSTUSEQS:
        if keepgoing:
            if len(seq_F) < len(lostuseq):
                break
            if seq_F[:len(lostuseq)] == lostuseq:
                seq_F = seq_F[len(lostuseq):]
            #if LOSTUSEQ[0] found, also look for LOSTUSEQ[1] etc.
            else:
                keepgoing = 0
    return seq_F


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument('afterBarcode', nargs='?', default=False)
    parser.add_argument('--SCA2', action='store_true', default=False)
    args = parser.parse_args()

    main(args.infile, args.afterBarcode, args.SCA2)
