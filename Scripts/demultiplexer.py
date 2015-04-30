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
#   To use, first hard-code the fixed sequence(s) below.
#
#######################################################################


import argparse
import time
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator


# Hard-code the fixed sequence(s) to look for after the barcode.
AFTER_BARCODE = ['ATCACGC',
                 'CGATGTC',
                 'TTAGGCC',
                 'TGACCAC',
                 'ACAGTGC',
                 'GCCAATC',
                 'CAGATCC',
                 'ACTTGAC',
                 'GATCAGC',
                 'TAGCTTC',
                 'GGCTACC',
                 'CTTGTAC',
                 'AGTCAAC',
                 'AGTTCCC',
                 'ATGTCAC',
                 'CCGTCCC',
                 'GTAGAGC',
                 'GTCCGCC',
                 'GTGAAAC',
                 'GTGGCCC',
                 'GTTTCGC',
                 'CGTACGC',
                 'GAGTGGC',
                 'GGTAGCC',
                 'ACTGATC',
                 'ATGAGCC',
                 'ATTCCTC',
                 'CAAAAGC',
                 'CAACTAC',
                 'CACCGGC']

BARCODE_LENGTH = 16
BARCODE_TRUNCATE = 0


def main(infile):
    forward_paired, reverse_paired, forward_unpaired = infile.split(',')
    readCount = 0
    start_time = time.time()
    i = 1
    for index in AFTER_BARCODE:
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
    args = parser.parse_args()

    main(args.infile)
