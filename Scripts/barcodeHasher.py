#!/usr/bin/env python

#############################################################################
#   barcodeHasher.py
#   2015 James A. Stapleton
#
#   This program sorts short reads into a dictionary on the basis of barcode
#       sequences. Output is a JSON-dumped file called unpairedDict.txt or
#       pairedDict.txt, depending on whether barcode pairing is used.
#
#   Arguments:
#       infile: comma-separated list of sequence files.
#
#       afterBarcode: one or more (comma-separated) sequences, one of
#               which is required to follow the barcode to confirm that
#               the read is not spurious
#
#   Options:
#
#       --SCA2: call when adapters with the SCA2 PCR primer
#               sequence are used
#               (SCA2: 5'-ACACGACGTGAACGATAGGAATTG-3')
#
#       --pairSeparateFile: call for two-tube barcode pairing
#
#       --pairSameFile: call for one-tube barcode pairing
#
#       --useFwdUnpaired: call to use forward reads whose reverse pairs
#                           were dropped by trimmomatic
#
#       --partial: print dictionaries every 1/NUMFRACTIONS of the data
#
#       --FLASH: merge overlapping forward and reverse reads with FLASH
#
#       --quality: add fastq quality line to the dictinary along with the
#                    sequence line, to allow error-checking by the assembler
#                       (e.g., SPAdes)
#
#############################################################################

import argparse
import time
import subprocess
import collections
import json
import itertools
import copy
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def main(infile, afterBarcode, SCA2, pairSeparateFile, pairSameFile,
         useFwdUnpaired, partial, FLASH, quality):

    if pairSeparateFile + pairSameFile > 1:
        print "Please choose only one barcode pairing option"
        return 1

    BARCODE_LENGTH = 16
    BARCODE_TRUNCATE = 1
    PAIR_THRESHOLD = 1
    ENDTRIM = 2
    NUMFRACTIONS = 10

    COMPLEMENT_DICT = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N'}
    AFTER_BARCODE_SEQS = []
    AFTER_BARCODE_RC_SEQS = []
    if afterBarcode:
        for ABseq in afterBarcode.split(','):
            AFTER_BARCODE_SEQS.append(ABseq)
            AFTER_BARCODE_RC_SEQS.append(''.join([COMPLEMENT_DICT[base]
                                         for base in ABseq])[::-1])
            print "After barcode:", ABseq
        ABcounts = [0] * len(AFTER_BARCODE_SEQS)
    else:
        ABcounts = [0]
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

    # take trimmed infiles, entered as file1.fq,file2.fq
    if pairSeparateFile and useFwdUnpaired:
        infile_F, infile_R, infile_unpaired, paired_F = infile.split(',')
    elif pairSameFile and useFwdUnpaired:
        infile_F, infile_R, infile_unpaired, trimmed = infile.split(',')
    elif pairSeparateFile:
        infile_F, infile_R, paired_F = infile.split(',')
    elif pairSameFile:
        infile_F, infile_R, trimmed = infile.split(',')
    elif useFwdUnpaired:
        infile_F, infile_R, infile_unpaired = infile.split(',')
    else:
        infile_F, infile_R = infile.split(',')

    #run FLASH to combine overlapping read pairs
    if FLASH:
        subprocess.call(["flash", "-M", "140", "-t", "1", infile_F, infile_R])

    # initilize stuff
    master_hash = collections.defaultdict(list)
    final_dict = collections.defaultdict(list)
    readCount = 0
    bad_after = 0
    passed = 0
    start_time = time.time()
    pairing = 0
    pair_dict = 0

# Pair barcodes
    if (pairSeparateFile + pairSameFile):
        pairing = 1
        print 'Pairing barcodes'

        # initilize stuff
        readCount = 0
        pair_dict = collections.defaultdict(lambda: collections.Counter())
        start_time = time.time()

        if pairSeparateFile:
            trimmed = paired_F

        with open(trimmed, 'rU') as merged:
            for title, seq, qual_F in FastqGeneralIterator(merged):
                # print time elapsed every 100000 reads processed
                readCount += 1
                if readCount % 100000 == 0:
                    print readCount
                    print time.time() - start_time, "seconds"
                    start_time = time.time()

                if len(seq) < BARCODE_LENGTH + BARCODE_LENGTH:
                    continue

                # pull out barcode1
                barcode1 = seq[BARCODE_TRUNCATE:BARCODE_LENGTH]
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
                seq_F, qual_F = trim_lost_U(seq_F, qual_F, LOSTUSEQS)
                for lostuseq_RC in LOSTUSEQ_RCS[::-1]:
                    lostarray = [lostuseq_RC]
                    seq_F, qual_F = trim_lost_U(seq_F, qual_F, lostarray)

                # check for RC of after_barcode_sequences
                if AFTER_BARCODE_SEQS:
                    test = after_barcode_seq(seq_F, AFTER_BARCODE_RC_SEQS)
                    if test == 1:
                        continue
                    else:
                        seq_F = seq_F[len(test):]

                # pull out barcode2
                barcode2rc = seq_F[: BARCODE_LENGTH - BARCODE_TRUNCATE]
                if len(barcode2rc) != BARCODE_LENGTH - BARCODE_TRUNCATE:
                    continue
                barcode2 = ''.join([COMPLEMENT_DICT[base] for 
                                    base in barcode2rc])[::-1]
                seq_F = seq_F[BARCODE_LENGTH:]

                pair_dict[barcode1][barcode2] += 1
                pair_dict[barcode2][barcode1] += 1

        #print pair_dict to file
        with open('barcodePairs.txt', 'w') as barcodePairs:
            print >> barcodePairs, json.dumps(pair_dict, sort_keys=True,
                                              indent=4, separators=(',', ': '))

    print 'Building the dictionary'
    readCount = 0

    # count the number of input reads to allow partial dictionaries
    numReads = 0
    if partial:
        if FLASH:
            with open("out.extendedFrags.fastq", 'rU') as merged:
                for line in merged:
                    numReads += 1
            with open("out.notCombined_1.fastq", 'rU') as unmerged_F:
                for line in unmerged_F:
                    numReads += 1
        else:
            with open(infile_F, 'rU') as unmerged_F:
                for line in unmerged_F:
                    numReads += 1
        if useFwdUnpaired:
            with open(infile_unpaired, 'rU') as unpaired_F:
                for line in unpaired_F:
                    numReads += 1
        numReads = numReads/4

    # iterate over merged reads from FLASH
    if FLASH:
        with open("out.extendedFrags.fastq", 'rU') as merged:
            for title, seq, qual in FastqGeneralIterator(merged):
                readCount, bad_after, ABcounts, passed, start_time = processReads(readCount, seq, qual, 0, 0, BARCODE_TRUNCATE,
                                                                                  BARCODE_LENGTH, ENDTRIM, bad_after, ABcounts,
                                                                                  AFTER_BARCODE_SEQS, AFTER_BARCODE_RC_SEQS,
                                                                                  LOSTUSEQS, passed, master_hash, pair_dict,
                                                                                  partial, numReads, NUMFRACTIONS, COMPLEMENT_DICT,
                                                                                  LOSTUSEQ_RCS, start_time,
                                                                                  PAIR_THRESHOLD, pairing, quality)

        # iterate over unmerged reads from FLASH
        with open("out.notCombined_1.fastq", 'rU') as unmerged_F:
            with open("out.notCombined_2.fastq", 'rU') as unmerged_R:
                f_iter = FastqGeneralIterator(unmerged_F)
                r_iter = FastqGeneralIterator(unmerged_R)
                for (title, seq, qual), (title_R, seq_R, qual_R) in itertools.izip(f_iter, r_iter):
                    readCount, bad_after, ABcounts, passed, start_time = processReads(readCount, seq, qual, seq_R, qual_R, BARCODE_TRUNCATE,
                                                                                      BARCODE_LENGTH, ENDTRIM, bad_after, ABcounts,
                                                                                      AFTER_BARCODE_SEQS, AFTER_BARCODE_RC_SEQS,
                                                                                      LOSTUSEQS, passed, master_hash, pair_dict,
                                                                                      partial, numReads, NUMFRACTIONS, COMPLEMENT_DICT,
                                                                                      LOSTUSEQ_RCS, start_time,
                                                                                      PAIR_THRESHOLD, pairing, quality)

    else:    # if not running FLASH, iterate over infiles
        with open(infile_F, 'rU') as unmerged_F:
            with open(infile_R, 'rU') as unmerged_R:
                f_iter = FastqGeneralIterator(unmerged_F)
                r_iter = FastqGeneralIterator(unmerged_R)
                for (title, seq, qual), (title_R, seq_R, qual_R) in itertools.izip(f_iter, r_iter):
                    readCount, bad_after, ABcounts, passed, start_time = processReads(readCount, seq, qual, seq_R, qual_R, BARCODE_TRUNCATE,
                                                                                      BARCODE_LENGTH, ENDTRIM, bad_after, ABcounts,
                                                                                      AFTER_BARCODE_SEQS, AFTER_BARCODE_RC_SEQS,
                                                                                      LOSTUSEQS, passed, master_hash, pair_dict,
                                                                                      partial, numReads, NUMFRACTIONS, COMPLEMENT_DICT,
                                                                                      LOSTUSEQ_RCS, start_time,
                                                                                      PAIR_THRESHOLD, pairing, quality)

    # iterate over forward_unpaired
    if useFwdUnpaired:
        print "Using forward unpaired reads"
        with open(infile_unpaired, 'rU') as unpaired_F:
            for title, seq, qual in FastqGeneralIterator(unpaired_F):
                readCount, bad_after, ABcounts, passed, start_time = processReads(readCount, seq, qual, 0, 0, BARCODE_TRUNCATE,
                                                                                  BARCODE_LENGTH, ENDTRIM, bad_after, ABcounts,
                                                                                  AFTER_BARCODE_SEQS, AFTER_BARCODE_RC_SEQS,
                                                                                  LOSTUSEQS, passed, master_hash, pair_dict,
                                                                                  partial, numReads, NUMFRACTIONS, COMPLEMENT_DICT,
                                                                                  LOSTUSEQ_RCS, start_time,
                                                                                  PAIR_THRESHOLD, pairing, quality)

    print str(readCount) + " total read pairs, " + str(passed) + " passed, " \
        + str(bad_after) + " non-compliant barcodes, AB counts: " + str(ABcounts)

    # print master_hash to file
    if (not pairSeparateFile) and (not pairSameFile):
        with open('unpairedDict.txt', 'w') as unpairedDict:
            print >> unpairedDict, json.dumps(master_hash, sort_keys=True,
                                              indent=4, separators=(',', ': '))
        del master_hash
        return 0

    else:
        # combine master_hash and pair_dict into final_dict
        with open('confirmedPairs.txt', 'w') as confirmed:
            final_dict = pairBarcodes(master_hash, pair_dict, final_dict,
                                      PAIR_THRESHOLD, confirmed)

        # print final_dict to file
        with open('pairedDict.txt', 'w') as pairedDict:
            print >> pairedDict, json.dumps(final_dict, sort_keys=True, indent=4,
                                            separators=(',', ': '))

    del master_hash
    del pair_dict
    del final_dict

    return 0


######## Function definitions ##############

def processReads(readCount, seq, qual_F, seq_R, qual_R, BARCODE_TRUNCATE,
                 BARCODE_LENGTH, ENDTRIM, bad_after, ABcounts,
                 AFTER_BARCODE_SEQS, AFTER_BARCODE_RC_SEQS, LOSTUSEQS,
                 passed, master_hash, pair_dict, partial, numReads,
                 NUMFRACTIONS, COMPLEMENT_DICT, LOSTUSEQ_RCS, start_time,
                 PAIR_THRESHOLD, pairing, quality):

    # Print time elapsed every 100000 reads processed
    readCount += 1
    if readCount % 100000 == 0:
        print readCount
        print time.time() - start_time, "seconds"
        start_time = time.time()

    # if testing partial fractions of the dataset,
    # print dicts every 1/NUMFRACTIONS of the data
    if partial:
        if not readCount % (numReads/NUMFRACTIONS):
            filenum = readCount/(numReads/NUMFRACTIONS)
            final_dict = collections.defaultdict(list)
            with open('unpairedDict'+str(filenum)+'.txt', 'w') as unpairedDict:
                print >> unpairedDict, json.dumps(master_hash, sort_keys=True,
                                                  indent=4, separators=(',', ': '))
            with open('confirmedPairs'+str(filenum)+'.txt', 'w') as confirmed:
                final_dict = pairBarcodes(master_hash, pair_dict, final_dict,
                                          PAIR_THRESHOLD, confirmed)
            if pairing:
                with open('pairedDict'+str(filenum)+'.txt', 'w') as pairedDict:
                    print >> pairedDict, json.dumps(final_dict, sort_keys=True,
                                                    indent=4, separators=(',', ': '))
                del final_dict
                final_dict = collections.defaultdict(list)

    barcode = seq[BARCODE_TRUNCATE:BARCODE_LENGTH]
    seq_F = seq[BARCODE_LENGTH:-ENDTRIM]
    qual_F = qual_F[BARCODE_LENGTH:-ENDTRIM]

    if seq_R and len(seq_R) > 30:
        barcode_RC = ''.join([COMPLEMENT_DICT[base] for base in barcode])[::-1]

        # look for RC of barcode at the end of a reverse read of a
        # fragment shorter than the read length
        if seq_R[-BARCODE_LENGTH : -BARCODE_TRUNCATE] == barcode_RC:
            seq_R = seq_R[ENDTRIM : -BARCODE_LENGTH]
            qual_R = qual_R[ENDTRIM : -BARCODE_LENGTH]
            # look for RC of the after-barcode sequences and lost-U sequences
            already_found_afterbarcodeseq = 0
            for afterbarcodeseqrc in AFTER_BARCODE_RC_SEQS:
                if already_found_afterbarcodeseq:
                    break
                test = 0
                test, seq_R, qual_R = reversetrim(seq_R, qual_R, afterbarcodeseqrc)
                if test:
                    for lostuseqrc in LOSTUSEQ_RCS:
                        if test:
                            test, seq_R, qual_R = reversetrim(seq_R, qual_R, lostuseqrc)
                        else:
                            already_found_afterbarcodeseq = 1
                            break
        else:
            seq_R = seq_R[ENDTRIM:-ENDTRIM]
            qual_R = qual_R[ENDTRIM:-ENDTRIM]

    test, seq_F, qual_F, seq_R, qual_R, bad_after, readCount, ABcounts = filters(seq_F, qual_F, seq_R, qual_R,
                                                                                 bad_after, readCount, ABcounts,
                                                                                 AFTER_BARCODE_SEQS, LOSTUSEQS)
    if test == 1:
        return readCount, bad_after, ABcounts, passed, start_time
    else:
        if seq_R and len(seq_R) > 30:
            seq_R, qual_R = sequence_checker(seq_R, qual_R, LOSTUSEQS)
            master_hash[barcode].append(seq_R)
            if quality:
                master_hash[barcode].append(qual_R)
        else:
            master_hash[barcode].append("")
            if quality:
                master_hash[barcode].append("")
        seq_F, qual_F = sequence_checker(seq_F, qual_F, LOSTUSEQS)
        master_hash[barcode].append(seq_F)
        if quality:
            master_hash[barcode].append(qual_F)
        passed += 1

        return readCount, bad_after, ABcounts, passed, start_time


def reversetrim(seq_R, qual_R, checkseq):
    '''trims sequences off front end of a reverse sequence read'''
    if seq_R[-len(checkseq):] == checkseq:
        seq_R = seq_R[:-len(checkseq)]
        qual_R = qual_R[:-len(checkseq)]
        return 1, seq_R, qual_R
    else:
        return 0, seq_R, qual_R


def filters(seq_F, qual_F, seq_R, qual_R, bad_after, readCount, ABcounts,
            AFTER_BARCODE_SEQS, LOSTUSEQS):
    "Collection of calls to the various filter functions"

    # read must have correct defined sequence after the barcode,
    # otherwise count and throw out
    if AFTER_BARCODE_SEQS:
        check = after_barcode_seq(seq_F, AFTER_BARCODE_SEQS)
        if check == 1:
            bad_after += 1
            return 1, seq_F, qual_F, seq_R, qual_R, bad_after, readCount, ABcounts
        else:
            i = 0
            for ABseq in AFTER_BARCODE_SEQS:
                if check == ABseq:
                    ABcounts[i] += 1
                    seq_F = seq_F[len(ABseq):]
                    qual_F = qual_F[len(ABseq):]
                    break
                i += 1

    # trim at N's
    seq_F, qual_F, seq_R, qual_R = Ntest(seq_F, qual_F, seq_R, qual_R)

    # look for lostU sequence after barcode, remove it
    seq_F, qual_F = trim_lost_U(seq_F, qual_F, LOSTUSEQS)

    # if everything is good,
    return 0, seq_F, qual_F, seq_R, qual_R, bad_after, readCount, ABcounts


def Ntest(seq_F, qual_F, seq_R, qual_R):
    "trim sequences with N's"
    Ntest = 0
    for i in xrange(0, len(seq_F)):
        if seq_F[i] == 'N':
            Ntest = 1
            break
    if Ntest == 1:
        seq_F = seq_F[0:i-1]
        qual_F = qual_F[0:i-1]

    if seq_R != 0:
        for i in xrange(1, len(seq_R)):
            if seq_R[i] == 'N':
                Ntest = 1
                break
        if Ntest == 1:
            seq_R = seq_R[0:i-1]
            qual_R = qual_R[0:i-1]
    return seq_F, qual_F, seq_R, qual_R


def after_barcode_seq(seq_F, AFTER_BARCODE_SEQS):
    "throw out sequences without the right sequence after the barcode"
    for ABseq in AFTER_BARCODE_SEQS:
        if seq_F[:len(ABseq)] == ABseq:
            return ABseq
    # if none of the allowed sequences are found,
    return 1


def sequence_checker(sequence, qual, LOSTUSEQS):
    "check the sequence for the LostU sequence or its RC and trim"
    # look for LOSTUSEQ or its RC with no mismatches. Trim.
    COMPLEMENT_DICT = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
    testseq = ''
    for lostuseq in LOSTUSEQS:
        testseq += lostuseq
        if len(testseq) > 7:
            testseq_RC = ''.join([COMPLEMENT_DICT[base] for base in testseq])[::-1]
            for i in xrange(len(sequence)-len(testseq)+1):
                if sequence[i:i+len(testseq)] == testseq or sequence[i:i+len(testseq)] == testseq_RC:
                    sequence = sequence[:i]
                    qual = qual[:i]
            break
    return sequence, qual


def trim_lost_U(seq_F, qual_F, LOSTUSEQS):
    """ test for lost U at the 3' end of the PCR primer sequence """
    keepgoing = 1
    for lostuseq in LOSTUSEQS:
        if keepgoing:
            if len(seq_F) < len(lostuseq):
                return seq_F, qual_F
            if seq_F[:len(lostuseq)] == lostuseq:
                seq_F = seq_F[len(lostuseq):]
                qual_F = qual_F[len(lostuseq):]
            #if LOSTUSEQ[0] found, also look for LOSTUSEQ[1] etc.
            else:
                keepgoing = 0
    return seq_F, qual_F


def pairBarcodes(master_hash, pair_dict, final_dict, PAIR_THRESHOLD, confirmed):
    numPairs = 0
    pair_dict_copy = copy.deepcopy(pair_dict)
    for barcode in master_hash:
        if barcode in pair_dict_copy:
             # sort the inner dict by values and take the 1st one = best_match
            sorted_pair_candidates = pair_dict_copy[barcode].most_common(1)
            if sorted_pair_candidates:
                best_match = sorted_pair_candidates[0][0]
                if best_match == 'X':
                    continue
                         # verify pair match by checking whether barcode
                         # is best_match's most frequent pair
                if (sorted_pair_candidates[0][1] > PAIR_THRESHOLD and
                        pair_dict_copy[best_match]):
                    # make sure best_match is f,
                    # i.e., hasn't been deleted in a prior cycle
                    cross_check = pair_dict_copy[best_match].most_common(1)
                    if (barcode == cross_check[0][0] and
                            cross_check[0][1] > PAIR_THRESHOLD and
                            master_hash.get(best_match) is not None):
                         # we have a verified pair!
                        numPairs += 1
                        confirmed.write(barcode + ' ' + best_match + '\n')
                        seq_list1 = master_hash.get(barcode)
                        seq_list2 = master_hash.get(best_match)
                        final_dict[barcode] = seq_list1 + seq_list2
                        # insert dummy as flag so when best_match is tried
                        # there is no duplication
                        pair_dict_copy[best_match]['X'] += 9999999
                    else:
                        final_dict[barcode] = master_hash[barcode]
                else:
                    final_dict[barcode] = master_hash[barcode]
            else:
                final_dict[barcode] = master_hash[barcode]
        else:
            final_dict[barcode] = master_hash[barcode]
    print "Found " + str(numPairs) + " barcode pairs"
    return final_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('afterBarcode', nargs='?', default=False)
    parser.add_argument('--SCA2', action='store_true', default=False)
    parser.add_argument('--pairSeparateFile', action='store_true', default=False)
    parser.add_argument('--pairSameFile', action='store_true', default=False)
    parser.add_argument('--useFwdUnpaired', action='store_true', default=False)
    parser.add_argument('--partial', action='store_true', default=False)
    parser.add_argument('--FLASH', action='store_true', default=False)
    parser.add_argument('--quality', action='store_true', default=False)
    args = parser.parse_args()

    main(args.infile, args.afterBarcode, args.SCA2, args.pairSeparateFile,
         args.pairSameFile, args.useFwdUnpaired, args.partial, args.FLASH, args.quality)
