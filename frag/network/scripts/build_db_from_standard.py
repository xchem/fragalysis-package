#!/usr/bin/env python
#
# Based on build.db.py, this module build the network from
# the Informatics Matters 'standard' (uncompressed) file representation.
#
# Alan Christie
# February 2019

import argparse
import os
import sys

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network, write_data


def main():
    """Read in a a standard file - then write out into a specified directory
    """
    parser = argparse.ArgumentParser(
        description="Convert un-compressed standard SMILES"
                    " for Astex Fragment network."
    )
    parser.add_argument("--input")
    parser.add_argument("--base_dir")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", dest="verbosity", action="store_const", const=1)
    group.add_argument("-vv", dest="verbosity", action="store_const", const=2)

    parser.set_defaults(verbosity=0)
    args = parser.parse_args()

    # Do we have an input and base directory?
    if not args.input:
        print('ERROR: Must specify an input')
        sys.exit(1)
    if not os.path.isfile(args.input):
        print('ERROR: input (%s) does not exist' % args.input)
        sys.exit(1)
    if not args.base_dir:
        print('ERROR: Must specify a base directory')
        sys.exit(1)
    if not os.path.isdir(args.base_dir):
        print('ERROR:input base directory (%s) does not exist' % args.base_dir)
        sys.exit(1)

    attrs = []
    standards = parser.parse_standard(args.input)
    for standard in standards:
        attrs.append(Attr(standard.smiles, ['EM', standard.compd_id]))

    # Build the network
    node_holder = NodeHolder()
    node_holder = build_network(attrs, node_holder,
                                args.base_dir,
                                args.verbosity)
    # Write the data out
    write_data(args.base_dir, node_holder, attrs)


if __name__ == "__main__":
    main()
