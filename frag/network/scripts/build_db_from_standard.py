#!/usr/bin/env python
#
# Based on build_db.py, this module builds the graph network from
# the Informatics Matters 'standard' (uncompressed) file representation.
#
# Alan Christie
# February 2019

import argparse
import os
import sys

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network, write_data
from frag.utils.parser import parse_standard_file


def main():
    """Read in a 'standard' file - then write out into a specified directory
    """
    parser = argparse.ArgumentParser(
        description="Convert un-compressed standard SMILES"
                    " for Astex Fragment network."
    )
    parser.add_argument("--input")
    parser.add_argument("--base_dir")
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise')
    parser.add_argument('--min-hac',
                        type=int, default=0,
                        help='Limit processing to molecules with at least this'
                             ' number of heady atoms')
    parser.add_argument('--max-hac',
                        type=int, default=0,
                        help='Limit processing to molecules with no more than'
                             ' this number of heavy atoms')
    parser.add_argument("--isomeric", dest="iso_flag", action="store_true")
    parser.add_argument("--non_isomeric", dest="iso_flag", action="store_false")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", dest="verbosity", action="store_const", const=1)
    group.add_argument("-vv", dest="verbosity", action="store_const", const=2)

    parser.set_defaults(verbosity=0)
    parser.set_defaults(iso_flag=True)
    args = parser.parse_args()

    # Do we have an input and base directory?
    if not args.input:
        print('ERROR: Must specify an input')
        sys.exit(1)
    if not os.path.isfile(args.input):
        print('ERROR: input (%s) does not exist' % args.input)
        sys.exit(2)
    if not args.base_dir:
        print('ERROR: Must specify a base directory')
        sys.exit(3)

    attrs = []
    standard_representations = parse_standard_file(args.input,
                                                   args.limit,
                                                   args.min_hac,
                                                   args.max_hac,
                                                   args.iso_flag)
    for standard_representation in standard_representations:
        attrs.append(Attr(standard_representation.smiles, []))

    if not os.path.isdir(args.base_dir):
        os.mkdir(args.base_dir)

    # Build the network
    node_holder = NodeHolder(iso_flag=args.iso_flag)
    node_holder = build_network(attrs, node_holder,
                                args.base_dir,
                                args.verbosity)
    # Write the data out
    write_data(args.base_dir, node_holder, attrs)


if __name__ == "__main__":
    main()
