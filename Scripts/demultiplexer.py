#!/usr/bin/env python

#######################################################################
#
#   demultiplexer.py
#   2015 James Stapleton
#
#   This script splits reads from a multiplexed sample prep
#   into separate FASTQ files on the basis of post-barcode indexes.
#
#   The algorithm looks for a fixed sequence(s) after a barcode.
#
#   Inputs:
#
#   infile: forward,reverse
#
#   afterBarcode: String of comma-separated index sequences to
#                    look for after the barcode
#
#   Options:
#
#   --BARCODE_LENGTH: Length of the barcode. Default 16.
#
#   --BARCODE_TRUNCATE: Ignore this many bases at the beginning of
#                       the barcode. Default 0.
#
#######################################################################


import argparse
import time
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def main(infile, afterBarcode, BARCODE_LENGTH, BARCODE_TRUNCATE):
    forward_paired, reverse_paired, forward_unpaired = infile.split(',')
    AFTER_BARCODE_SEQS = []
    for ABseq in afterBarcode.split(','):
        AFTER_BARCODE_SEQS.append(ABseq)
    readCount = 0
    start_time = time.time()
    i = 1
    for index in AFTER_BARCODE_SEQS:
        with open(forward_paired, 'rU') as fwd_paired:
            with open(reverse_paired, 'rU') as rev_paired:
                f_iter = FastqGeneralIterator(fwd_paired)
                r_iter = FastqGeneralIterator(rev_paired)
                with open("forward_paired_" + str(i) + ".fq", 'a') as outfileF:
                    with open("reverse_paired_" + str(i) + ".fq", 'a') as outfileR:
                        for (title, seq, qual), (title_R, seq_R, qual_R) in itertools.izip(f_iter, r_iter):
                            # print time elapsed every 1000000 reads processed
                            readCount += 1
                            if readCount % 1000000 == 0:
                                print readCount
                            # check for AFTER sequence
                            after_length = len(index)
                            if seq[BARCODE_LENGTH : BARCODE_LENGTH +
                                   after_length] == index:
                                outfileF.write('@' + title + '\n' +
                                               seq + '\n+\n' + qual + '\n')
                                outfileR.write('@'+ title_R + '\n' +
                                               seq_R + '\n+\n' + qual_R + '\n')

        with open(forward_unpaired, 'rU') as fwd_unpaired:
            with open("forward_unpaired_" + str(i) + ".fq", 'a') as outfile:
                for title, seq, qual in FastqGeneralIterator(fwd_unpaired):
                    # print time elapsed every 100000 reads processed
                    readCount += 1
                    if readCount % 1000000 == 0:
                        print readCount
                    # check for AFTER sequence
                    if seq[BARCODE_LENGTH : BARCODE_LENGTH + after_length] == index:
                        outfile.write('@' + title + '\n' + seq + '\n+\n' + qual + '\n')

        print "Finished with index", str(i)
        print time.time() - start_time, "seconds"
        start_time = time.time()
        readCount = 0
        i += 1

    print str(readCount) + " total reads"

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument('afterBarcode', nargs='?', default=['ATCACGC,CGATGTC,TTAGGCC,TGACCAC,ACAGTGC,GCCAATC,CAGATCC,ACTTGAC,GATCAGC,TAGCTTC,GGCTACC,CTTGTAC,AGTCAAC,AGTTCCC,ATGTCAC,CCGTCCC,GTAGAGC,GTCCGCC,GTGAAAC,GTGGCCC,GTTTCGC,CGTACGC,GAGTGGC,GGTAGCC,ACTGATC,ATGAGCC,ATTCCTC,CAAAAGC,CAACTAC,CACCGGC'],
        help='one or more (comma-separated) sequences, one of which is required to follow the barcode to confirm that the read is not spurious')
    parser.add_argument('--BARCODE_LENGTH', action="store", dest="BARCODE_LENGTH", type=int, default=16,
            help='length of the barcode, default 16.')
    parser.add_argument('--BARCODE_TRUNCATE', action="store", dest="BARCODE_TRUNCATE", type=int, default=0,
            help='Ignore this many bases at the beginning of the barcode, default 0.')
    args = parser.parse_args()

    main(args.infile, args.afterBarcode, args.BARCODE_LENGTH, args.BARCODE_TRUNCATE)
