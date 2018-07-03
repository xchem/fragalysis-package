import argparse

from rdkit import Chem
from tqdm import tqdm
from frag.utils.parser import get_file, parse_mols

if __name__ == "__main__":

    # Read in a SD or SMILES file - then write out into a specified directory
    parser = argparse.ArgumentParser(
        description="Convert a SMILES or SDFile to input for Astex Fragment network."
    )
    parser.add_argument("--input")
    parser.add_argument("--input_format", default="smi")
    parser.add_argument("--output")
    parser.add_argument("--output_format", default="smi")
    parser.add_argument("--chunk_size", default=None)
    parser.add_argument("--num_chunks", default=None)
    args = parser.parse_args()

    if args.chunk_size:
        chunk_size = int(args.chunk_size)
        num_chunks = None
    elif args.num_chunks:
        num_chunks = int(args.num_chunks)
        chunk_size = None

    mols = parse_mols(args.input, args.input_format)

    if num_chunks:
        out_files = [
            get_file(args.output, args.output_format, x) for x in range(num_chunks)
        ]
    counter = 0
    file_counter = 0
    for x in tqdm(mols):
        # Ignore None mols
        if x is None:
            continue
        if num_chunks:
            out_file = out_files[counter % num_chunks]
        else:
            if counter % chunk_size == 0:
                out_file = get_file(args.output, args.output_format, file_counter)
                file_counter += 1
        out_file.write(x)
        counter += 1
