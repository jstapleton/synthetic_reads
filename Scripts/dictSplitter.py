#!/usr/bin/env python

##################################################################
#
#   dictSplitter.py
#   2015 James Stapleton
#
#   Splits a read dictionary into multiple smaller ones
#       that can be assembled in parallel
#
##################################################################

import argparse


def main(infile, num_smaller_dicts):
    with open(infile, 'r') as data:
        barcode_count = 0
        for line in data:
            if line[-2] == '[':
                barcode_count += 1

    with open(infile, 'r') as data:
        barcodes_per_file = barcode_count/int(num_smaller_dicts)
        this_file_count = 0
        seq_list = []
        i = 1
        outfile = open("splitDict" + str(i) + ".txt", "w")
        for line in data:
            if line[-2] == '}':
                break
            elif line[0] == '{':
                seq_list.append(line)
            elif line[-2] == '[':
                seq_list.append(line)
            elif line[-3] == '"':
                seq_list.append(line)
            elif line[-2] == '"':
                seq_list.append(line)
            elif line[-3] == ']' or line[-2] == ']':
                seq_list.append(line)
                if this_file_count > barcodes_per_file:
                    outfile.write("}\n")
                    outfile.close()
                    i += 1
                    print i
                    outfile = open("splitDict" + str(i) + ".txt", "w")
                    this_file_count = 0
                    outfile.write("{\n")
                for line in seq_list:
                    outfile.write(line)
                this_file_count += 1
                seq_list = []

    outfile.write("}\n")
    outfile.close()

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('num_smaller_dicts')
    args = parser.parse_args()

    main(args.infile, args.num_smaller_dicts)
