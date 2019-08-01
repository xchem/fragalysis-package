#!/usr/bin/env python

"""deduplicate_chemspace_bb_compounds.py

Given a ChemSpace file this module writes out the file with
unique compound IDs.

We create a '.dedup.gz' file that contains a 1st-line
'header' formed from the _OUTPUT_COLUMNS list.

Alan Christie
August 2019
"""

import argparse
import glob
import gzip
import logging
import sys

# Configure basic logging
logger = logging.getLogger('chemspace')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The minimum number of columns in the input files and
# and a map of expected column names indexed by (0-based) column number.
#
# The 'standardised' files contain at least 2 columns...
#
# smiles    0
# id        1

expected_min_num_cols = 2
smiles_col = 0
compound_col = 1
expected_input_cols = {compound_col: 'id',
                       smiles_col: 'smiles'}

# All the vendor compound IDs
vendor_compounds = set()
# A map of duplicate compounds and the number of duplicates.
# The index uses the vendor's original ID value, not our prefixed value.
duplicate_suffix = '-duplicate-'
vendor_duplicates = {}


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


def deduplicate_vendor_compounds(output_file, file_name):
    """Process the given file and deduplicate the compound IDs.

    :param output_file: The deduplicated output file
    :param file_name: The (compressed) file to process
    :returns: The number of items processed
    """
    global vendor_compounds
    global vendor_duplicates
    global duplicate_suffix

    logger.info('Deduplicating %s...', file_name)

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

            fields = line.split('\t')
            if len(fields) <= 1:
                continue

            osmiles = fields[smiles_col].strip()
            vendor_id = fields[compound_col].strip()

            # Add the compound (expected to be unique)
            # to our set of 'all compounds'.
            deduplicated_vendor_id = vendor_id
            if deduplicated_vendor_id in vendor_compounds:
                # Get the number of duplicates (default of 1)
                # using the vendor's original ID as a key
                duplicate_count = vendor_duplicates.setdefault(vendor_id, 1)
                deduplicated_vendor_id = '{}{}{}'.format(vendor_id,
                                                         duplicate_suffix,
                                                         duplicate_count)
                # Increment for any further duplicates
                vendor_duplicates[vendor_id] = duplicate_count + 1
            else:
                vendor_compounds.add(vendor_id)

            # Write the deduplicated data

            output = [osmiles,
                      deduplicated_vendor_id]

            output_file.write('\t'.join(output) + '\n')

    return num_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Deduplicator (ChemSpace/BB)')
    parser.add_argument('vendor_dir',
                        help='The ChemSpace vendor directory,'
                             ' containing tab-delimited the ".gz" files'
                             ' to be processed.')
    parser.add_argument('vendor_prefix',
                        help='The ChemSpace vendor file prefix,'
                             ' i.e. "Jul2019". Only files with this prefix'
                             ' in the vendor directory will be processed')
    parser.add_argument('output',
                        help='The output basename')

    args = parser.parse_args()

    # Output is either s fixed name in an output directory
    # or a prefixed filename (without an output directory)
    output_filename = '{}.dedupe.gz'.format(args.output)

    # Before we open the output file
    # get a list of all the input files (the prefix may be the same)
    # so we don't want our file in the list of files to be processed)
    data_files = glob.glob('{}/{}*.gz'.format(args.vendor_dir,
                                              args.vendor_prefix))

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    logger.info('Writing %s...', output_filename)
    num_processed = 0
    with gzip.open(output_filename, 'wt') as output_gzip_file:

        # Write the header...
        output_gzip_file.write('\t'.join(['SMILES', 'ID']) + '\n')

        # Process all the Vendor files...
        for data_file in data_files:
            num_processed += deduplicate_vendor_compounds(output_gzip_file,
                                                          data_file)

    # Summary
    for vendor_duplicate in vendor_duplicates:
        logger.info('Duplicate compound: {} x{}'.
                    format(vendor_duplicate,
                           vendor_duplicates[vendor_duplicate]))
    logger.info('{:,} compounds with duplicates'.format(len(vendor_duplicates)))
