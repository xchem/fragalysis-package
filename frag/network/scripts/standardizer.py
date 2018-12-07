#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import glob, gzip, sys, argparse
from frag.utils.rdkit_utils import standardize


uncharger = rdMolStandardize.Uncharger()

def process(smiles):

    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None, 0
        frags = Chem.GetMolFrags(m, asMols=True)
        mol = standardize(m)
        return mol, len(frags)
    except ValueError as ve:
        sys.stderr.write("Failed to process molecule\n")
        sys.stderr.write(ve.message)

    return None, 0


def log(sep, vals):
    sys.stderr.write(sep.join(vals) + '\n')

def next_output_file_name(basename, chunksize, counter):
    if chunksize:
        filename = basename + "_" + str(counter) + ".smi"
    else:
        filename = basename + ".smi"
    sys.stderr.write("Opening output file " + filename + "\n")
    return filename

parser = argparse.ArgumentParser(description='Fragment network standardize')
parser.add_argument('-o', '--output', help='Base name for output', default='output')
parser.add_argument('-l', '--limit', help='Only process this many from each file', type=int)
parser.add_argument('-c', '--chunk-size', help='Generate multiple output files with this many records', type=int)
# to use tab as a delimiter set command line option as -d $'\t'
parser.add_argument('--delimiter-input', help='record delimter for input data', default=' ')
parser.add_argument('--delimiter-output', help='record delimter for output data', default=' ')
parser.add_argument('--min-hac', help='Minimum HAC', type=int, default=0)
parser.add_argument('--max-hac', help='Maximum HAC', type=int)
parser.add_argument('--smiles-column', help='The molecule SMILES column (first column is 0)', default=0, type=int)
parser.add_argument('--id-column', help='The molecule ID column (first column is 0)', required=True, type=int)
parser.add_argument('--id-prefix', help='The molecule ID prefix to insert', required=True)
parser.add_argument('inputs', nargs='+')
args = parser.parse_args()

outfilename = args.output
delimiter_input = args.delimiter_input
delimiter_output = args.delimiter_output
chunk_size = args.chunk_size # can be None
if chunk_size:
    log(' ', ["Using chunk size of", str(chunk_size)])

inputFileCount = 0
outputFileCount = 0
outfile = None

for file in args.inputs:
    prefix_delimiter = ':'
    countMols = 0
    excluded = 0
    included = 0
    errors = 0

    inputFileCount += 1
    log(' ', ['Opening input file', str(inputFileCount), file])

    with gzip.open(file, 'rb') as f:
        linecount = 0
        for line in f:
            linecount += 1
            # skip the header line
            if linecount == 1:
                continue
            if args.limit and linecount > args.limit + 1:
                break

            values = line.strip().split(delimiter_input)
            #sys.stderr.write("Found " + str(len(values))  + " values")
            if len(values) == 0:
                continue
            countMols += 1

            if outfile is None or (chunk_size > 0 and (countMols % chunk_size) == 0):
                if outfile is not None:
                    outfile.close()
                outputFileCount += 1
                outfile = open(next_output_file_name(outfilename, chunk_size, outputFileCount), "w")
                outfile.write("SSMILES ID OSMILES\n")

            osmiles = values[0]
            m, numFrags = process(osmiles)
            if m is not None:
                ssmiles = Chem.MolToSmiles(m)
                hac = m.GetNumHeavyAtoms()
                #print("HAC: " + str(hac))
                vals = [ssmiles,
                        '{}{}{}'.format(args.id_prefix,
                                        prefix_delimiter,
                                        values[args.id_column]),
                        osmiles]
                if hac >= args.min_hac and (args.max_hac is None or hac <= args.max_hac):
                    included += 1
                    outfile.write(delimiter_output.join(vals) + "\n")
                    outfile.flush()

                else:
                    excluded += 1
            else:
                errors += 1

    log('\n', [
        str(included) + ' included',
        str(excluded) + ' excluded',
        str(errors) + ' errors\n'])

outfile.close()