#!/usr/bin/env python
#
# Filters a standard file based on skip/limit and min/max HAC.
#
# Alan Christie
# May 2019

import argparse
import os
import sys

from frag.std_utils.parser import filter_standard_file


def main():
    """Read in a 'standard' file - printing out a smaller version
    """
    parser = argparse.ArgumentParser(
        description="Filter compressed standard SMILES"
                    " for Astex Fragment network."
    )
    parser.add_argument('--input',
                        help="The input, a compressed standard file.")
    parser.add_argument('--output',
                        help="The filtered output file, normally ending in"
                             " 'tab.gz' and written as a compressed file")
    parser.add_argument('--reject',
                        help="The rejected output file containing molecules"
                             "considered but not selected, normally ending in"
                             " 'tab.gz' and written as a compressed file")
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise')
    parser.add_argument('-s', '--skip',
                        type=int, default=0,
                        help='Number of molecules to skip molecules'
                             ' in the input file')
    parser.add_argument('--min-hac',
                        type=int, default=0,
                        help='Limit processing to molecules with at least this'
                             ' number of heavy atoms')
    parser.add_argument('--max-hac',
                        type=int, default=0,
                        help='Limit processing to molecules with no more than'
                             ' this number of heavy atoms')

    args = parser.parse_args()

    # Do we have an input and base directory?
    if not args.input:
        print('ERROR: Must specify an input')
        sys.exit(1)
    if not os.path.isfile(args.input):
        print('ERROR: input (%s) does not exist' % args.input)
        sys.exit(2)

    filter_standard_file(args.input,
                         args.limit,
                         args.skip,
                         args.min_hac,
                         args.max_hac,
                         args.output,
                         args.reject)


if __name__ == "__main__":
    main()
