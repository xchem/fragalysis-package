#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import glob, gzip, sys, argparse


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


def standardize(mol):
    mol = rdMolStandardize.Cleanup(mol)
    mol = fragment(mol);
    mol = uncharger.uncharge(mol)
    return mol


def fragment(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)

    if len(frags) == 1:
        return mol
    else:
        # TODO - handle ties
        biggest_index = -1
        i = 0

        biggest_count = 0
        for frag in frags:
            hac = frag.GetNumHeavyAtoms()
            if hac > biggest_count:
                biggest_count = hac
                biggest_mol = frag
                biggest_index = i
            i+=1

        return biggest_mol


def log(sep, vals):
    sys.stderr.write(sep.join(vals) + '\n')

parser = argparse.ArgumentParser(description='Fragment network standardize')
parser.add_argument('-o', '--output', help='Base name for output', default='output')
parser.add_argument('-l', '--limit', help='Only process this many from each file', type=int)
parser.add_argument('--min-hac', help='Minimum HAC', type=int, default=0)
parser.add_argument('--max-hac', help='Maximum HAC', type=int)
parser.add_argument('--id-column', help='The molecule ID column (first column is 1)', type=int)
parser.add_argument('--id-prefix', help='The molecule ID prefix to insert')
parser.add_argument('inputs', nargs='+')
args = parser.parse_args()

if args.id_column < 1:
    sys.stderr.write("ID column must be 1 or higher\n")
    sys.exit(1)

outfilename = args.output

prefix_delimiter = ':'
countMols = 0
countFiles = 0
excluded = 0
included = 0
errors = 0

for file in args.inputs:
    countFiles += 1
    log(' ', [str(countFiles), file])
    outfile = open(outfilename + "_" + str(countFiles) + ".smi", "w")
    outfile.write("SSMILES OSMILES ID\n")
    with gzip.open(file, 'rb') as f:
        linecount = 0
        for line in f:
            linecount += 1
            # skip the header line
            if linecount == 1:
                continue
            if args.limit and linecount > args.limit + 1:
                break
            countMols += 1
            values = line.split(" ")
            osmiles = values[0]
            m, numFrags = process(osmiles)
            if m is not None:
                ssmiles = Chem.MolToSmiles(m)
                hac = m.GetNumHeavyAtoms()
                #print("HAC: " + str(hac))
                vals = [ssmiles, osmiles, '{}{}{}'.format(args.id_prefix,
                                                          prefix_delimiter,
                                                          values[args.id_column])]
                if hac >= args.min_hac and (args.max_hac is None or hac <= args.max_hac):
                    included += 1
                    outfile.write(" ".join(vals) + "\n")
                    outfile.flush()

                else:
                    excluded += 1
            else:
                errors += 1

    outfile.close()

    log('\n', [
        str(included) + ' included',
        str(excluded) + ' excluded',
        str(errors) + ' errors\n'])
