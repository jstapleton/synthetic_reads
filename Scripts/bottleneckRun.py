#!/usr/bin/env python

#####################################################################
# This script is for use when running a MiSeq run to
# determine the complexity of a sample
# before preparing a library for synthetic long read assembly.
#
# The algorithm looks for a fixed sequence after a barcode,
# then counts the number of
# times in the file each barcode precedes that sequence.
#
# To use, first hard-code the fixed sequence below.
#
#
#####################################################################


import argparse
import time
import collections
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from scipy import stats

####### Hard-code the fixed sequence(s) to look for after the barcode.

AFTER_BARCODE = ['ATCACG',
                 'CGATGT',
                 'TTAGGC',
                 'TGACCA',
                 'ACAGTG',
                 'GCCAAT',
                 'CAGATC',
                 'ACTTGA',
                 'GATCAG',
                 'TAGCTT',
                 'GGCTAC',
                 'CTTGTA',
                 'AGTCAA',
                 'AGTTCC',
                 'ATGTCA',
                 'CCGTCC',
                 'GTAGAG',
                 'GTCCGC',
                 'GTGAAA',
                 'GTGGCC',
                 'GTTTCG',
                 'CGTACG',
                 'GAGTGG',
                 'GGTAGC',
                 'ACTGAT',
                 'ATGAGC',
                 'ATTCCT',
                 'CAAAAG',
                 'CAACTA',
                 'CACCGG']

BARCODE_LENGTH = 16
BARCODE_TRUNCATE = 0  # disregard the first positions in a barcode
THRESHOLD = 2


def main(infile):

    # initilize stuff
    readCount = 0
    start_time = time.time()
    unaccounted = 0

    # create a data structure to count barcode
    # occurrences for each AFTER sequence:
    counter_list = []
    for sequence in AFTER_BARCODE:
        counter_list.append(collections.Counter())

    with open(infile, 'rU') as fwdReads:
        # iterate over forward reads looking for paired barcodes
        for title, seq, qual in FastqGeneralIterator(fwdReads):
            # print time elapsed every 100000 reads processed
            readCount += 1
            if readCount % 100000 == 0:
                print readCount
                print time.time() - start_time, "seconds"
                start_time = time.time()

            # check for AFTER sequences
            i = 0
            for sequence in AFTER_BARCODE:
                after_length = len(sequence)
                if seq[BARCODE_LENGTH : BARCODE_LENGTH + after_length] == sequence:
                    barcode = seq[BARCODE_TRUNCATE:BARCODE_LENGTH]
                    counter_list[i][barcode] += 1
     #               break
                if i == len(AFTER_BARCODE) - 1:
                    unaccounted += 1
                i += 1

    for counter in counter_list:
        above_threshold = 0
        for barcode in counter:
            if counter[barcode] >= THRESHOLD:
                above_threshold += 1
        print above_threshold

#    # return a histogram list of number of barcodes seen 1x, 2x, etc for each index
#    for counter in counter_list:
#        print sum(counter.values())  # first print total number of reads for barcode
#        for num in xrange(1,10):
#            at_num = 0
#            for barcode in counter:
#                if counter[barcode] == num:
#                   at_num += 1
#            print at_num
#        print '\n'

#    i = 0
#    for counter in counter_list:
#        with open('counter'+str(i), 'w') as outfile:
#            print >> outfile, json.dumps(counter, sort_keys=True,
#                                         indent=4, separators=(',', ': '))
#        i += 1

#    print str(readCount) + " total reads"
#    print str(unaccounted) + " reads unaccounted for"

#    with open('histogram3.txt', 'w') as outfile:
#        for barcode in counter_list[2]:
#            outfile.write(str(counter_list[2][barcode]) + '\n')

    # create a histogram list of number of barcodes seen 1x, 2x, etc for each index
    for counter in counter_list:
        total_reads = sum(counter.values())     # total number of reads for barcode
        histogram = [0]
        for num in range(1, 10):
            at_num = 0
            for barcode in counter:
                if counter[barcode] == num:
                    at_num += 1
            histogram.append(at_num)

        # error correction: ~8% of barcodes have errors,
        # remove this many from total count and from number of barcodes seen 1x
        histogram[1] = histogram[1] - total_reads * 0.08
        total_reads = total_reads * 0.92

        # normalize histogram to total number of barcodes
        total_barcodes = 0
        for i in range(1, len(histogram)):
            total_barcodes += histogram[i]
        for i in range(1, len(histogram)):
            if total_barcodes:
                histogram[i] = histogram[i] / float(total_barcodes)

        # find the poisson mu that minimizes the square difference to the histogram
        best = 1000
        mu = 0
        for x in range(1, 40):
            test_mu = x * 0.1
            sum_diff_sq = 0
            for i in range(2, 5):
                # remove possibility of x = 0 since we can't see that in data
                poisson_value = (stats.distributions.poisson.pmf(i, test_mu)
                                 / (1 - stats.distributions.poisson.pmf(0, test_mu)))
                diff = histogram[i] - poisson_value
                diff_sq = diff * diff
                sum_diff_sq += diff_sq

            if sum_diff_sq < best:
                best = sum_diff_sq
                mu = test_mu

        molecules = total_reads / mu
#        print histogram
#        print total_reads
#        print mu
        print molecules
   #     print '\n'

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    args = parser.parse_args()

    main(args.infile)
