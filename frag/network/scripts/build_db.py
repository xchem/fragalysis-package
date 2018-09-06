import argparse
import os

from rdkit import Chem

from frag.network.models import NodeHolder, Attr
from frag.utils.network_utils import build_network, write_data
from frag.utils.parser import parse_mols

from tqdm import tqdm


def main():
    # Read in a SD or SMILES file - then write out into a specified directory
    parser = argparse.ArgumentParser(
        description="Convert a SMILES or SDFile to input for Astex Fragment network."
    )
    parser.add_argument("--input")
    parser.add_argument("--input_format", default="smi")
    parser.add_argument("--base_dir")

    args = parser.parse_args()
    attrs = []
    id = 0
    mols = parse_mols(args.input, args.input_format)
    for x in tqdm(mols):
        if x is None:
            continue
        attr = Attr(
            Chem.CanonSmiles(Chem.MolToSmiles(x, isomericSmiles=True)),
            ["EM", x.GetProp("_Name")],
        )
        attrs.append(attr)
        id += 1
    if not os.path.isdir(args.base_dir):
        os.mkdir(args.base_dir)
    # Build the network
    node_holder = NodeHolder()
    node_holder = build_network(attrs, node_holder)
    # Write the data out
    write_data(args.base_dir, node_holder, attrs)


if __name__ == "__main__":
    main()
