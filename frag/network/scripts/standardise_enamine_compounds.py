#!/usr/bin/env python

"""standardise_real_compounds.py

Processes Enamine Real vendor compound files, and generates a 'standard'
tab-separated output.

We create a 'real-standardised-compounds.tab' file that contains a 1st-line
'header' formed from the _OUTPUT_COLUMNS list.

Alan Christie
January 2019
"""

import argparse
import glob
import gzip
import logging
import os
import sys

from rdkit import RDLogger

from standardise_molport_compounds import standardise

# Configure basic logging
logger = logging.getLogger('real')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The columns in our output file.
_OUTPUT_COLUMNS = ['OSMILES',
                   'ISO_SMILES',
                   'NONISO_SMILES',
                   'CMPD_ID']

# The minimum number of columns in the input files and
# and a map of expected column names indexed by (0-based) column number.
#
# The 'standardised' files contain at least 3 columns...
#
# smiles    0
# idnumber  1

expected_min_num_cols = 2
smiles_col = 0
compound_col = 1
expected_input_cols = {compound_col: 'idnumber',
                       smiles_col: 'smiles'}

# Prefix for output files
output_filename_prefix = 'real'
# The output file.
# Which will be gzipped.
output_filename = output_filename_prefix + '-standardised-compounds.tab'

# The compound identifier prefix
# the vendor uses in the the compound files...
supplier_prefix = 'Z'
# The prefix we use in our fragment file
real_prefix = 'REAL:'

# All the vendor compound IDs
vendor_compounds = set()

# Various diagnostic counts
num_vendor_mols = 0
num_vendor_molecule_failures = 0

# The line rate at which the process writes updates to stdout.
# Every 1 million?
report_rate = 1000000


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


def standardise_vendor_compounds(output_file, file_name):
    """Process the given file and standardise the vendor
    information, writing it as tab-separated fields to the output.

    As we load the Vendor compounds we 'standardise' the SMILES and
    determine whether they represent an isomer or not.

    :param output_file: The tab-separated standardised output file
    :param file_name: The (compressed) file to process
    """
    global vendor_compounds
    global num_vendor_mols
    global num_vendor_molecule_failures

    logger.info('Standardising %s...', file_name)

    line_num = 0
    with gzip.open(file_name, 'rt') as gzip_file:

        # Check first line (a space-delimited header).
        # This is a basic sanity-check to make sure the important column
        # names are what we expect.

        hdr = gzip_file.readline()
        field_names = hdr.split()
        # Expected minimum number of columns...
        if len(field_names) < expected_min_num_cols:
            error('expected at least {} columns found {}'.
                  format(expected_input_cols, len(field_names)))
        # Check salient columns (ignoring case)...
        for col_num in expected_input_cols:
            actual_name = field_names[col_num].strip().lower()
            if actual_name != expected_input_cols[col_num]:
                error('expected "{}" in column {} found "{}"'.
                      format(expected_input_cols[col_num],
                             col_num,
                             actual_name))

        # Columns look right...

        for line in gzip_file:

            line_num += 1
            fields = line.split()
            if len(fields) <= 1:
                continue

            if line_num % report_rate == 0:
                logger.info(' ...at compound {:,}'.format(line_num))

            osmiles = fields[smiles_col]
            compound_id = real_prefix + fields[compound_col]

            # Add the compound (expected to be unique)
            # to our set of 'all compounds'.
            if compound_id in vendor_compounds:
                error('Duplicate compound ID ({})'.format(compound_id))
            vendor_compounds.add(compound_id)

            # Standardise and update global maps...
            # And try and handle and report any catastrophic errors
            # from dependent modules/functions.

            std, iso, noniso = standardise(osmiles)
            if not std:
                num_vendor_molecule_failures += 1
                continue
            num_vendor_mols += 1

            # Write the standardised data

            output = [osmiles,
                      iso,
                      noniso,
                      compound_id]

            output_file.write('\t'.join(output) + '\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Standardiser (MolPort)')
    parser.add_argument('vendor_dir',
                        help='The Enamine vendor directory,'
                             ' containing the ".gz" files to be processed.')
    parser.add_argument('vendor_prefix',
                        help='The Enamine vendor file prefix,'
                             ' i.e. "June2018". Only files with this prefix'
                             ' in the vendor directory will be processed')
    parser.add_argument('output',
                        help='The output directory')

    args = parser.parse_args()

    # Create the output directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if not os.path.isdir(args.output):
        error('output ({}) is not a directory'.format(args.output))

    # -------
    # Stage 1 - Process Vendor Files
    # -------

    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    output_filename = os.path.join(args.output, '{}.gz'.format(output_filename))
    logger.info('Writing %s...', output_filename)
    with gzip.open(output_filename, 'wt') as output_gzip_file:

        # Write the header...
        output_gzip_file.write('\t'.join(_OUTPUT_COLUMNS) + '\n')

        # Process all the Vendor files...
        real_files = glob.glob('{}/{}*.gz'.format(args.vendor_dir,
                                                     args.vendor_prefix))
        for real_file in real_files:
            standardise_vendor_compounds(output_gzip_file, real_file)

    # Summary
    logger.info('{:,} vendor molecules'.format(num_vendor_mols))
    logger.info('{:,} vendor molecule failures'.format(num_vendor_molecule_failures))
