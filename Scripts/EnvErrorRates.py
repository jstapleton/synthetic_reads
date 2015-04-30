#!/usr/bin/env python

######################################################################
#
#   EnvErrorRates.py
#   2015 James Stapleton
#
#   Parses error statistics grepped from EMBOSS water,
#       prints error rate
#
######################################################################

import argparse


def main(errorRates):

    with open(errorRates, 'r') as contigs:
        contiglist = contigs.readlines()

    errorlist = []
    for line in contiglist:
      #  if line[1:11] == 'Identities':
        #if line[1:11] == 'Identity: ':
            #slashline = line.split()[2]
            #errorlist.append(slashline)

    #for slash in errorlist:
        #num = slash.split('/')[0]
        #denom = slash.split('/')[1]
        if line == '\n':
            print '0'
            continue
        num = line.split('/')[0]
        denom = line.split('/')[1]
       # print denom

        print float(num)/float(denom)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('errorRates')
    args = parser.parse_args()

    main(args.errorRates)
