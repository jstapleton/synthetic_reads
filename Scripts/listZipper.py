#!/usr/bin/env python

##################################################################
#
#   listZipper.py
#   2015 James Stapleton
#
#   Compares two lists of numbers and writes a new list
#       with the larger number for each pair
#
##################################################################


import argparse
from itertools import izip


def main(infile1, infile2):
    with open(infile1, 'r') as fwd, open(infile2, 'r') as RC:
        for x, y in izip(fwd, RC):
            if float(x) > float(y):
                print float(x)
            else:
                print float(y)
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile1')
    parser.add_argument('infile2')
    args = parser.parse_args()

    main(args.infile1, args.infile2)
