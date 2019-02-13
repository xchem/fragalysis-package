#!/usr/bin/env python
# coding=utf-8

"""standardise_molport_compounds.py

Processes MolPort vendor compound files, expected to contain pricing
information and generates a 'standard' tab-separated output.
We create a 'molport-standardised-compounds.tab' file that contains a 1st-line
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

from frag.utils import standardise_utils

# The columns in our output file.
_OUTPUT_COLUMNS = standardise_utils.STANDARD_COLUMNS + \
                  ['PRICERANGE_1MG',
                   'PRICERANGE_5MG',
                   'PRICERANGE_50MG',
                   'BEST_LEAD_TIME']

# Configure basic logging
logger = logging.getLogger('molport')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The minimum number of columns in the input data and
# a map of expected column names indexed by column number.
#
# The molecule data is spread over a number of `txt.gz` files
# (i.e. files like `iis_smiles-000-000-000--000-499-999.txt.gz`)
# in a common directory where the files have the following header
# names and (0-based) positions:
#
# SMILES                0
# SMILES_CANONICAL      1
# MOLPORTID             2
# STANDARD_INCHI        3
# INCHIKEY              4
# PRICERANGE_1MG        5
# PRICERANGE_5MG        6
# PRICERANGE_50MG       7
# BEST_LEAD_TIME        8

expected_min_num_cols = 9
smiles_col = 0
compound_col = 2
cost_col = {1: 5, 5: 6, 50: 7}
blt_col = 8
expected_input_cols = {smiles_col: 'SMILES',
                       compound_col: 'MOLPORTID',
                       cost_col[1]: 'PRICERANGE_1MG',
                       cost_col[5]: 'PRICERANGE_5MG',
                       cost_col[50]: 'PRICERANGE_50MG',
                       blt_col: 'BEST_LEAD_TIME'}

# The output file.
# Which will be gzipped.
output_filename = 'standardised-compounds.tab'

# The compound identifier prefix
# the vendor uses in the the compound files...
supplier_prefix = 'MolPort-'
# The prefix we use in our fragment file
# and the prefix we use for our copy of the
molport_prefix = 'MOLPORT:'

# All the vendor compound IDs
vendor_compounds = set()

# Various diagnostic counts
num_vendor_mols = 0
num_vendor_molecule_failures = 0

# The line rate at which the process writes updates to stdout.
report_rate = 250000


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: %s', msg)
    sys.exit(1)


def standardise_vendor_compounds(output_file, file_name, limit):
    """Process the given file and standardise the vendor (and pricing)
    information, writing it as tab-separated fields to the output.

    As we load the Vendor compounds we 'standardise' the SMILES and
    determine whether they represent an isomer or not.

    :param output_file: The tab-separated standardised output file
    :param file_name: The (compressed) file to process
    :param limit: Limit processing to this number of values (or all if 0)
    :returns: The number of items processed
    """

    global vendor_compounds
    global num_vendor_mols
    global num_vendor_molecule_failures

    logger.info('Standardising %s...', file_name)

    line_num = 0
    num_processed = 0
    with gzip.open(file_name, 'rt') as gzip_file:

        # Check first line (a space-delimited header).
        # This is a basic sanity-check to make sure the important column
        # names are what we expect.

        hdr = gzip_file.readline()
        field_names = hdr.split('\t')
        # Expected minimum number of columns...
        if len(field_names) < expected_min_num_cols:
            error('expected at least {} columns found {}'.
                  format(expected_input_cols, len(field_names)))
        # Check salient columns...
        for col_num in expected_input_cols:
            actual_name = field_names[col_num].strip().upper()
            if actual_name != expected_input_cols[col_num]:
                error('expected "{}" in column {} found "{}"'.
                      format(expected_input_cols[col_num],
                             col_num,
                             actual_name))

        # Columns look right...

        for line in gzip_file:

            line_num += 1
            fields = line.split('\t')
            if len(fields) <= 1:
                continue

            if line_num % report_rate == 0:
                logger.info(' ...at compound {:,}'.format(line_num))

            osmiles = fields[smiles_col]
            compound_id = molport_prefix + fields[compound_col].split(supplier_prefix)[1]

            # Add the compound (expected to be unique)
            # to our set of 'all compounds'.
            if compound_id in vendor_compounds:
                error('Duplicate compound ID ({})'.format(compound_id))
            vendor_compounds.add(compound_id)

            # Standardise and update global maps...
            # And try and handle and report any catastrophic errors
            # from dependent modules/functions.

            std_info = standardise_utils.standardise(osmiles)
            if not std_info.std:
                num_vendor_molecule_failures += 1
                continue
            num_vendor_mols += 1

            # Write the standardised data

            output = [osmiles,
                      std_info.iso,
                      std_info.noniso,
                      std_info.hac,
                      compound_id,
                      fields[cost_col[1]].strip(),
                      fields[cost_col[5]].strip(),
                      fields[cost_col[50]].strip(),
                      fields[blt_col].strip()]

            output_file.write('\t'.join(output) + '\n')

            # Enough?
            num_processed += 1
            if limit and num_processed >= limit:
                break

    return num_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Standardiser (MolPort)')
    parser.add_argument('vendor_dir',
                        help='The MolPort vendor directory,'
                             ' containing the ".gz" files to be processed.')
    parser.add_argument('vendor_prefix',
                        help='The MolPort vendor file prefix,'
                             ' i.e. "iis_smiles". Only files with this prefix'
                             ' in the vendor directory will be processed')
    parser.add_argument('output',
                        help='The output directory')
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise')

    args = parser.parse_args()

    # Create the output directory
    if os.path.exists(args.output):
        logger.error('Output exists')
        sys.exit(1)
    os.mkdir(args.output)
    os.chmod(args.output, 0o777)

    # -------
    # Stage 1 - Process Vendor Files
    # -------

    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Report any limiting...?
    if args.limit:
        logger.warning('Limiting processing to first {:,} molecules'.format(args.limit))

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    output_filename = os.path.join(args.output, '{}.gz'.format(output_filename))
    logger.info('Writing %s...', output_filename)
    num_processed = 0
    with gzip.open(output_filename, 'wt') as output_gzip_file:

        # Write the header...
        output_gzip_file.write('\t'.join(_OUTPUT_COLUMNS) + '\n')

        # Process all the Vendor files...
        molport_files = glob.glob('{}/{}*.gz'.format(args.vendor_dir,
                                                     args.vendor_prefix))
        for molport_file in molport_files:
            num_processed += standardise_vendor_compounds(output_gzip_file,
                                                          molport_file,
                                                          args.limit)
            if args.limit and num_processed >= args.limit:
                break

    # Summary
    logger.info('{:,} vendor molecules'.format(num_vendor_mols))
    logger.info('{:,} vendor molecule failures'.format(num_vendor_molecule_failures))
