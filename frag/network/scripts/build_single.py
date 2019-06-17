#!/usr/bin/env python
#
# Build nodes and edges for a single smiles


import argparse
import sys

from rdkit import Chem

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network
from frag.utils.rdkit_utils import standardize


def main():
    # Read in a SD or SMILES file - then write out into a specified directory
    parser = argparse.ArgumentParser(
        description="Convert a SMILES to nodes, edges and attributes"
    )
    parser.add_argument("--smiles")
    parser.add_argument("--id")
    parser.add_argument("--standardize", action="store_true")
    parser.add_argument("--isomeric", dest="iso_flag", action="store_true")
    parser.add_argument("--non_isomeric", dest="iso_flag", action="store_false")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", dest="verbosity", action="store_const", const=1)
    group.add_argument("-vv", dest="verbosity", action="store_const", const=2)

    parser.set_defaults(verbosity=0)
    parser.set_defaults(iso_flag=True)

    args = parser.parse_args()

    # Do we have an input and base directory?
    if not args.smiles:
        print('ERROR: Must specify a SMILES')
        sys.exit(1)

    attrs = []
    print("Original SMILES: " + args.smiles)
    mol = Chem.MolFromSmiles(args.smiles)
    if args.standardize:
        mol = standardize(mol)
        print("Standardized SMILES: " + Chem.MolToSmiles(mol))
    smiles = Chem.CanonSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
    print("Canonical SMILES: " + smiles)

    id = args.id
    if id is None:
        id = "smiles1"
    attr = Attr(smiles, ["EM", id])
    attrs.append(attr)
    # Build the network
    node_holder = NodeHolder(iso_flag=args.iso_flag)
    max_frags = 0
    node_holder = build_network(attrs, node_holder,
                                max_frags, smiles, args.verbosity)
    # Write the data out
    for node in node_holder.node_list:
        print(str(node))

    for edge in node_holder.get_edges():
        print(str(edge))

    for attr in attrs:
        print(str(attr))

    print("Number of nodes: " + str(len(node_holder.node_list)))
    print("Number of edges: " + str(len(node_holder.get_edges())))


if __name__ == "__main__":
    main()
