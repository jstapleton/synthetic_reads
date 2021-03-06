#!/usr/bin/env python

#############################################################################
#   dictReader.py
#   2015 James A. Stapleton
#
#   This program parses JSON-dumped barcode-hashed short-read dictionaries
#       created by barcodeHasher.py and sends files to SPAdes or Velvet
#       for assembly.
#
#   Options:
#       --makeHistogram: counts the number of reads in each  barcode-defined
#                           group, returns a file histogram.txt with one
#                           count per line that can be used to make a
#                           histogram and estimate the sample complexity.
#                           The constant MIN_NUMBER_OF_SEQUENCES defines the
#                           cutoff below which a bin is ignored.
#
#       --runVelvet: creates a fasta file from the reads in each barcode-
#                       defined group and sends it to Velvet for assembly.
#
#       --diginorm: performs digital normalization on each read group before
#                      sending it for assembly.
#
#       --runSpades: creates fasta or fastq files from the reads in each
#                       barcode-defined group and sends them to SPAdes for
#                       assembly.
#
#       --runTruSpades: creates fastq files from the reads in each
#                          barcode-defined group and sends them to truSPAdes
#                          for assembly.
#
#       --HPCC: calls SPAdes with one thread as required by some HPCCs.
#
#       --quality: creates fastq files rather than fasta files and uses
#                   SPAdes's built-in error correction mode.
#
#       --TRUNCATED_BARCODE_LENGTH: BARCODE_LENGTH - BARCODE_TRUNCATE from
#                                       the barcodeHasher.py run, default 16
#
#       --MIN_NUMBER_OF_SEQUENCES: minimum number of reads in a bin to count
#                                      in --makeHistogram
#
#       --threads: number of threads to use when running SPAdes, default 8,
#                  --HPCC sets to 1
#
#       --runCelera: calls the Celera assembler
#
#       --runTadpole: calls the Tadpole assember from bbmap
#
#############################################################################

from __future__ import division
from itertools import izip
import argparse
import subprocess
import os


def main(infile, makeHistogram, runVelvet, diginorm, runSpades, runTruSpades, HPCC,
        quality, TRUNCATED_BARCODE_LENGTH, MIN_NUMBER_OF_SEQUENCES, threads, runTadpole, runCelera):

    KMER_LENGTH = 99 # for diginorm
    MIN_CONTIG_LENGTH = 350 # for Velvet

    if not runVelvet and not makeHistogram and not runSpades and not runTruSpades and not runTadpole and not runCelera:
        print "Not doing anything!"
        return 0

    if makeHistogram:
        #histogram = collections.Counter()
        histoOut = open("histogram.txt", 'w')
        goodReads = 0
        totalReads = 0

    if (runVelvet + runSpades + runTruSpades):
        if os.path.exists('./contig_list.txt'):
            subprocess.call(['rm', 'contig_list.txt'])

    if runTruSpades:
        quality = 1
        if not os.path.exists('./truspades_input'):
            os.makedirs('./truspades_input')

    if runTadpole or runCelera:
        quality = 1

    with open(infile, 'r') as data:
        # parse JSON-dumped dictionary (uniqueDict.txt, unpairedDict.txt, etc)
        # Consider one barcode-defined group at a time
        seq_list = []
        for line in data:
            if line[-2] == '}':     # end of the file
                break
            elif line[0] == '{':    # start of the file
                continue
            # start of a barcode, example:     "AAAAAACGTTATGCAG": [
            elif line[-2] == '[':
                barcode = line[5:5+TRUNCATED_BARCODE_LENGTH]
            # start of a sequence, example:
            # "GGAAACTATACTAAAACTTGCTAAAAGCCATGATAAACTGAT",
            elif line[-3] == '"':
                sequence = line[9:-3]
                seq_list.append(sequence)
            # start of the last sequence in a group
            # no comma at the end, so -2 instead of -3
            elif line[-2] == '"':
                sequence = line[9:-2]
                seq_list.append(sequence)
            # end of a barcode-defined group, example:     ],
            elif line[-3] == ']' or line[-2] == ']':
                complete_list = seq_list
                seq_list = []

                if makeHistogram:
                    # this is the old way of printing a
                    # histogram of barcode abundance:
                    # returns pairs like 100, 20
                     # but this is hard to visualize
    #                 numBarcodes = len(complete_list)
    #                 histogram[numBarcodes] += 1
                    # new way: just print the number
                    # of sequences in the barcode
                    #  group to a file, plot a histogram with Python
                   # if len(complete_list) > 10 and len(complete_list) < 10000:
                    if quality:
                        if len(complete_list) > 2 * MIN_NUMBER_OF_SEQUENCES:
                            # divide by two here to get number of PE read pairs
                            print >> histoOut, len(complete_list)/4
                        totalReads += len(complete_list)/2
                        if len(complete_list) > 1000:
                            goodReads += 500
                        elif len(complete_list) > 200:
                            goodReads += len(complete_list)/2
                    else:
                        if len(complete_list) > MIN_NUMBER_OF_SEQUENCES:
                            # divide by two here to get number of PE read pairs
                            print >> histoOut, len(complete_list)/2
                        totalReads += len(complete_list)
                        if len(complete_list) > 500:
                            goodReads += 500
                        elif len(complete_list) > 100:
                            goodReads += len(complete_list)

                # write sequences from each items[barcode] to a fasta file,
                # send to assembler, overwrite fasta with next barcode sequnces
                # only run assembler if there are enough reads
                # for good coverage
                if (((not quality) and
                    (len(complete_list) > MIN_NUMBER_OF_SEQUENCES))
                    or (quality and
                        (len(complete_list) > 2 * MIN_NUMBER_OF_SEQUENCES))):

                    if runVelvet:
                            # open an output file that will be
                        # overwritten each time through
                        # the loop with a fasta list of sequences
                        # to feed to Velvet
                        with open('fasta.txt', 'w') as fasta:
                            i = 1
                            for seq in complete_list:
                                print >> fasta, ('>' + str(i))
                                if seq == '':
                                    fasta.write('A\n')
                                else:
                                    fasta.write(seq+'\n')
                                i += 1

                    # send fasta to Velveth
                    # Diginorm
                        if diginorm:
                            subprocess.call(["normalize-by-median.py",
                                             "-k", "20", "-C", "20", "-N", "4",
                                             "-x", "5e8", "fasta.txt"])
                            subprocess.call(["velveth", "seqData",
                                             str(KMER_LENGTH),
                                             "fasta.txt.keep"])
                        else:
                            subprocess.call(["velveth", "seqData",
                                             str(KMER_LENGTH), "-shortPaired",
                                             "fasta.txt"])

                        # get back files from Velveth
                        # send these files to Velvetg
                        subprocess.call(["velvetg", "seqData", "-exp_cov",
                                         "auto", "-cov_cutoff", 'auto',
                                         "-min_contig_lgth",
                                         str(MIN_CONTIG_LENGTH)])

                        # append contigs.fa to a growing file of contigs
                        with open("seqData/contigs.fa", "r") as fin:
                            data = fin.read()
                        with open("contig_list.txt", "a") as fout:
                            fout.write('>' + "Barcode: " + barcode + "\n")
                            fout.write(data + "\n")
                        with open("seqData/stats.txt", "r") as statsin:
                            stats = statsin.read()
                        with open("stats_list.txt", "a") as statsout:
                            statsout.write('>'+"Barcode: " + barcode + "\n")
                            statsout.write(stats + "\n")

                    if runSpades or runTadpole or runCelera:
                        # open an output file that will be
                        # overwritten each time through
                        # the loop with a fasta list of
                        # sequences to feed to Spades

                        if quality:
                            with open('left.fq', 'w') as left,\
                                    open('right.fq', 'w') as right,\
                                    open('unpaired.fq', 'w') as unpaired:
                                i = 1
                                for seq1, qual1, seq2, qual2 in\
                                        grouper(4, complete_list):
                                    if (seq1 == "") and (seq2 == ""):
                                        continue
                                    elif seq1 == "":
                                        unpaired.write('@Seq_ID' +
                                                       str(i) + '\n')
                                        unpaired.write(seq2 + '\n')
                                        unpaired.write('+\n')
                                        unpaired.write(qual2 + '\n')
                                    elif seq2 == "":
                                        unpaired.write('@Seq_ID' +
                                                       str(i) + '\n')
                                        unpaired.write(seq1 + '\n')
                                        unpaired.write('+\n')
                                        unpaired.write(qual1 + '\n')
                                    else:
                                        left.write('@Seq_ID' + str(i) + '\n')
                                        left.write(seq1 + '\n')
                                        left.write('+\n')
                                        left.write(qual1 + '\n')
                                        right.write('@Seq_ID' + str(i) +
                                                    '\n')
                                        right.write(seq2 + '\n')
                                        right.write('+\n')
                                        right.write(qual2 + '\n')
                                    i += 1

                        else:
                            with open('fasta.fa', 'w') as fasta:
                                i = 1
                                for seq in complete_list:
                                    print >> fasta, ('>' + str(i))
                                    if seq == "":
                                        print >> fasta, "A"
                                    else:
                                        print >> fasta, seq
                                    i += 1

                    if runSpades:
                    # send fasta to Spades
                        if HPCC:
                            threads = 1

                        if quality:
                            if diginorm:
                                subprocess.call(["normalize-by-median.py",
                                                 "-k", "21", "-C", "20",
                                                 "-N", "4", "-x", "5e8",
                                                 "-p", "-s",
                                                 "normC20k20.kh",
                                                 "fastq.fq"])
                                subprocess.call(["filter-abund.py", "-V",
                                                 "normC20k20.kh",
                                                 "fastq.fq.keep"])
                                subprocess.call(["mv", "fastq.fq.keep",
                                                 "fastq.fq"])
                                subprocess.call(["spades.py", "-k",
                                                 "21,33,55,77,99,127",
                                                 "-t", str(threads), "--careful",
                                                 "--sc",
                                                 "--pe1-12", "fastq.fq",
                                                 "-o", "spades_output"])
                            else:
                                subprocess.call(["spades.py", "-k",
                                                 "21,33,55,77,99,127",
                                                 "-t", str(threads),
                                                 "--careful", "--sc",
                                                 "--pe1-1", "left.fq",
                                                 "--pe1-2", "right.fq",
                                                 "--pe1-s", "unpaired.fq",
                                                 "-o", "spades_output",
                                                 "--disable-gzip-output"])
                        else:
                            if diginorm:
                                subprocess.call(["normalize-by-median.py",
                                                 "-k", "21", "-C", "20",
                                                 "-N", "4", "-x", "5e8",
                                                 "-p", "-s",
                                                 "normC20k20.kh",
                                                 "fasta.fa"])
                                subprocess.call(["filter-abund.py", "-V",
                                                 "normC20k20.kh",
                                                 "fast.fa.keep"])
                                subprocess.call(["mv", "fasta.fa.keep",
                                                 "fasta.fa"])
                                subprocess.call(["spades.py", "-k",
                                                 "21,33,55,77,99,127",
                                                 "-t", str(threads),
                                                 "--careful", "--sc",
                                                 "--pe1-12", "fasta.fa",
                                                 "-o", "spades_output"])
                            else:
                                subprocess.call(["spades.py", "-k",
                                                 "21,33,55,77,99,127",
                                                 "-t", str(threads), "--careful",
                                                 "--only-assembler",
                                                 "--sc",
                                                 "--pe1-12", "fasta.fa",
                                                 "-o", "spades_output"])

                        # append contigs.fasta to a growing file of contigs
                        if os.path.exists("./spades_output/contigs.fasta"):
                            with open("./spades_output/contigs.fasta",
                                      "r") as fin:
                                data = fin.read()
                            with open('contig_list.txt', 'a') as fout:
                                fout.write('>' + "Barcode: " + barcode + "\n")
                                fout.write(data + "\n")

                    if runTruSpades:
                        # open an output file that will be
                        # overwritten each time through
                        # the loop with a fastq list of
                        # sequences to feed to truSpades

                        with open('./truspades_input/reads_L1_R1.fastq', 'w') as left,\
                                open('./truspades_input/reads_L1_R2.fastq', 'w') as right:
                                i = 1
                                for seq1, qual1, seq2, qual2 in\
                                        grouper(4, complete_list):
                                    if (seq1 == "") or (seq2 == ""):
                                        continue
                                    left.write('@Seq_ID' + str(i) + '\n')
                                    left.write(seq1 + '\n')
                                    left.write('+\n')
                                    left.write(qual1 + '\n')
                                    right.write('@Seq_ID' + str(i) +
                                                '\n')
                                    right.write(seq2 + '\n')
                                    right.write('+\n')
                                    right.write(qual2 + '\n')
                                    i += 1

                                  # write dataset file
                        with open('dataset_file.txt', 'w') as dataset_file:
                            current_path = os.path.dirname(os.path.realpath('dataset_file.txt'))
                            dataset_file.write(barcode + ' ' + current_path + '/truspades_input/reads_L1_R1.fastq ' \
                                    + current_path + '/truspades_input/reads_L1_R2.fastq')

                    # run truSpades
                        subprocess.call(["truspades.py",
                                         "--dataset", "dataset_file.txt",
                                         "-t", str(threads),
                                         "-o", "truspades_output"])

                        # append contigs.fasta to a growing file of contigs
                        if os.path.exists("./truspades_output/TSLR.fasta"):
                            with open("./truspades_output/TSLR.fasta",
                                      "r") as fin:
                                data = fin.read()
                            with open('contig_list.txt', 'a') as fout:
                                fout.write('>' + "Barcode: " + barcode + "\n")
                                fout.write(data + "\n")

                    if runTadpole:
                        subprocess.call(["tadpole.sh",
                                          "in=left.fq",
                                          "in2=right.fq",
                                          "overwrite=true",
                                          "oute1=correct1.fq",
                                          "oute2=correct2.fq",
                                          "mode=extend",
                                          "el=50", "er=50",
                                          "k=62", "ecc=t"])
                        subprocess.call(["tadpole.sh",
                                          "in=correct1.fq",
                                          "in2=correct2.fq",
                                          "out=tadpole.fa",
                                          "overwrite=true",
                                          "k=90",
                                          "rinse=t", "shave=t",
                                          "mincontig=250",
                                          "mincov=3", "bm1=8"])

                         # append contigs.fasta to a growing file of contigs
                        if os.path.exists("tadpole.fa"):
                            with open("tadpole.fa", "r") as fin:
                                data = fin.read()
                            with open('contig_list_tadpole.txt', 'a') as fout:
                                fout.write('>' + "Barcode: " + barcode + "\n")
                                fout.write(data + "\n")

                    if runCelera:
                        with open('leftright.frg', 'w') as outfile:
                            subprocess.call(["fastqToCA", "-insertsize", "250", "25",
                                             "-libraryname", "celera", "-technology",
                                             "illumina", "-mates", "left.fq,right.fq",
                                             "-nonrandom"], stdout=outfile)
                        subprocess.call(["runCA", "-d", "celera", "-p", "celera", "leftright.frg"])

                        if os.path.exists("celera/9-terminator/celera.utg.fasta"):
                            with open("celera/9-terminator/celera.utg.fasta", "r") as fin:
                                data = fin.read()
                            with open('contig_list_celera.txt', 'a') as fout:
                                fout.write('>' + "Barcode: " + barcode + "\n")
                                fout.write(data + "\n")
                            subprocess.call(["rm", "-r", "celera", "leftright.frg"])

    if makeHistogram:
#         for numBarcodes in histogram:
#             print >> histoOut, (str(numBarcodes) + ','
#                           + str(histogram[numBarcodes]))
        histoOut.close()
#         del histogram
        print "efficiency score: ", goodReads/totalReads
        print "total number of reads: ", totalReads

    return 0


def grouper(n, iterable):
    "s -> (s0,s1,...sn-1), (sn,sn+1,...s2n-1), (s2n,s2n+1,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)

    # http://stackoverflow.com/questions/4356329/
    # creating-a-python-dictionary-from-a-line-of-text/4356415#4356415

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--makeHistogram', action='store_true', default=False,
            help='counts the number of reads in each barcode-defined group,\
            returns a file histogram.txt with one count per line that can be used to\
            make a histogram and estimate the sample complexity. The constant\
            MIN_NUMBER_OF_SEQUENCES defines the cutoff below which a bin is ignored.')
    parser.add_argument('--runVelvet', action='store_true', default=False,
            help='creates a fasta file from the reads in each barcode-defined group and sends it to Velvet for assembly.')
    parser.add_argument('--diginorm', action='store_true', default=False,
            help='performs digital normalization on each read group before sending it for assembly.')
    parser.add_argument('--runSpades', action='store_true', default=False,
            help='creates fasta or fastq files from the reads in each barcode-defined group and sends them to SPAdes for assembly.')
    parser.add_argument('--runTruSpades', action='store_true', default=False,
            help='creates fastq files from the reads in each barcode-defined group and sends them to truSPAdes for assembly.')
    parser.add_argument('--HPCC', action='store_true', default=False,
            help='calls SPAdes with one thread as required by some HPCCs.')
    parser.add_argument('--quality', action='store_true', default=False,
            help='creates fastq files rather than fasta files and uses SPAdes built-in error correction mode.')
    parser.add_argument('--TRUNCATED_BARCODE_LENGTH', action='store', dest="TRUNCATED_BARCODE_LENGTH", type=int, default=16,
            help='BARCODE_LENGTH - BARCODE_TRUNCATE from the barcodeHasher.py run, default 16')
    parser.add_argument('--MIN_NUMBER_OF_SEQUENCES', action='store', dest="MIN_NUMBER_OF_SEQUENCES", type=int, default=100,
            help='minimum number of reads in a bin to count in --makeHistogram')
    parser.add_argument('--threads', action='store', dest="threads", type=int, default=8,
            help='number of threads to use when running SPAdes, default 8, --HPCC sets to 1')
    parser.add_argument('--runTadpole', action='store_true', default=False,
            help='calls the Tadpole assembler from bbmap')
    parser.add_argument('--runCelera', action='store_true', default=False,
            help='calls the Celera assembler')
    args = parser.parse_args()

    main(args.infile, args.makeHistogram, args.runVelvet,
         args.diginorm, args.runSpades, args.runTruSpades,
         args.HPCC, args.quality, args.TRUNCATED_BARCODE_LENGTH,
         args.MIN_NUMBER_OF_SEQUENCES, args.threads, args.runTadpole, args.runTruSpades)
