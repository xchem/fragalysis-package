#!/usr/bin/env python

"""process_excluded_to_standard.py

Given an excluded file (from a prior build) this module
fabricates a standard file. It simply replicates the
SMILES and sets HAC to 0 and compound id to something unique for th row.

We create a 'standardised-compounds.tab' file that contains a 1st-line
'header' formed from the _OUTPUT_COLUMNS list.

Alan Christie
July 2019
"""

import argparse
import gzip
import logging
import os
import sys

from frag.std_utils import parser

# Configure basic logging
logger = logging.getLogger('process')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The columns in our output file.
# In this file we don't add any of our own.
_OUTPUT_COLUMNS = parser.STANDARD_COLUMNS

# The output file.
# Which will be gzipped.
output_filename = 'standardised-compounds.tab'

# The prefix we use in our fragment file
compound_prefix = 'fabricated-'


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


def standardise_excluded_compounds(output_file, file_name):
    """Process the given file and standardise the excluded content.
    The excluded file contains a SMILES string in the 2nd column
    (comma-separated).

    :param output_file: The tab-separated standardised output file
    :param file_name: The (compressed) file to process
    :returns: The number of items processed
    """
    line_num = 0
    lines_processed = 0
    with gzip.open(file_name, 'rt') as gzip_file:
        for line in gzip_file:

            line_num += 1
            fields = line.split(',')
            if len(fields) <= 1:
                continue

            smiles = fields[1].strip()

            # Write the standardised data
            output = [smiles,
                      smiles,
                      smiles,
                      '0',
                      '%s%s' % (compound_prefix, line_num)]
            output_file.write('\t'.join(output) + '\n')
            lines_processed += 1

    return lines_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Exclusion Standardiser')
    parser.add_argument('exclusion_file',
                        help='The file of exclusions.'
                             ' A ".gz" of excluded SMILES strings.')
    parser.add_argument('output',
                        help='The output directory, where'
                             ' the standard file will be written')

    args = parser.parse_args()

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    output_filename = os.path.join(args.output,
                                   '{}.gz'.format(output_filename))
    logger.info('Writing %s...', output_filename)
    num_processed = 0
    with gzip.open(output_filename, 'wt') as output_gzip_file:

        # Write the header...
        output_gzip_file.write('\t'.join(_OUTPUT_COLUMNS) + '\n')

        # Process the excluded file...
        num_processed = standardise_excluded_compounds(output_gzip_file,
                                                       args.exclusion_file)
    logger.info('Done (%s)', num_processed)
