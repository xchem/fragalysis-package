#!/usr/bin/env python

"""process_enamine_compounds.py

Processes standardised Enamine vendor files, not expected to contain
pricing information.

Two new files are generated and the original nodes file augmented with a
"V_E" label.

Note:   This module does expect `colate_all` to have been used on the original
        graph files to produce normalised supplier identities in the node file
        this module uses.

The purpose of this module is to create "Vendor" Compound nodes
and relationships to augment the DLS fragment database.
Every fragment line that has an Enamine identifier in the original data set
is labelled and a relationship created between it and the Vendor's compound(s).

Some vendor compound nodes may not exist in the original data set.

The files generated (in a named output directory) are:

-   "enamine-compound-nodes.csv.gz"
    containing all the nodes for the vendor compounds.

-   "enamine-molecule-compound_edges.csv.gz"
    containing the relationships between the original node entries and
    the "Vendor" nodes. There is a relationship for every Enamine
    compound that was found in the earlier processing.

The module augments the original nodes by adding the label
"V_E" for all MolPort compounds that have been found
to the augmented copy of the original node file that it creates.

If the original nodes file is "nodes.csv" the augmented copy
(in the named output directory) will be called
"enamine-augmented-nodes.csv.gz".

Alan Christie
November 2018
"""

import argparse
from collections import namedtuple
import glob
import gzip
import logging
import os
import re
import sys

# Configure basic logging
logger = logging.getLogger('enamine')
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
# The 'standardised' files contain at least 3 columns...
#
# smiles    0
# idnumber  1

expected_min_num_cols = 2
smiles_col = 0
compound_col = 1
expected_input_cols = {compound_col: 'idnumber',
                       smiles_col: 'smiles'}

# The Vendor Compound node has...
# a compound id (A 'Z' number)
# a SMILES string
CompoundNode = namedtuple('CompoundNode', 'c s')
# The unique set of vendor compound nodes.
compound_nodes = set()
# The vendor compound IDs that have a node.
# The index is the compound ID, the value is the VendorNode
compound_map = {}

# Prefix for output files
output_filename_prefix = 'enamine'
# The namespaces of the various indices
smiles_namespace = 'F2'
compound_namespace = 'VE'

# Regular expression to find
# Enamine compound IDs in the original nodes file.
# These can have a 'Z' or 'PV-' prefix.
enamine_re = re.compile(r'REAL:((?:Z|PV\-)\d+)')

# Various diagnostic counts
num_nodes = 0
num_nodes_augmented = 0
num_compound_relationships = 0


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


def extract_vendor(gzip_filename):
    """Process the given file and extract vendor information.
    """

    global compound_nodes

    logger.info('Processing {}...'.format(gzip_filename))

    num_lines = 0
    with gzip.open(gzip_filename, 'rt') as gzip_file:

        # Check first line (a space-delimited header).
        # This is a basic sanity-check to make sure the important column
        # names are what we expect.

        hdr = gzip_file.readline()
        field_names = hdr.split()
        # Expected minimum number of columns...
        if len(field_names) < expected_min_num_cols:
            error('expected at least {} columns found {}'.
                  format(expected_input_cols, len(field_names)))
        # Check salient columns...
        for col_num in expected_input_cols:
            if field_names[col_num].strip() != expected_input_cols[col_num]:
                error('expected "{}" in column {} found "{}"'.
                      format(expected_input_cols[col_num],
                             col_num,
                             field_names[col_num]))

        # OK - looks like the column names are right.
        # let's load the data...

        for line in gzip_file:

            num_lines += 1
            fields = line.split()

            smiles = fields[smiles_col]
            compound_id = fields[compound_col]

            # The compound ID must be unique
            if compound_id in compound_map:
                error('Duplicate compound "{}"'.format(compound_id))

            compound_node = CompoundNode(compound_id, smiles)
            compound_nodes.add(compound_node)
            # Put in a map, indexed by compound ID for fast lookup later
            compound_map[compound_id] = compound_node


def write_compound_nodes(directory, compounds):
    """Writes the CompoundNodes to a node file, including a header.

    :param directory: The output directory
    :param compounds: The set of compound namedtuples
    """

    filename = os.path.join(directory,
                            '{}-compound-nodes.csv.gz'.
                            format(output_filename_prefix))
    logger.info('Writing {}...'.format(filename))

    with gzip.open(filename, 'wb') as gzip_file:
        gzip_file.write('cmpd_id:ID({}),'
                        'smiles,'
                        ':LABEL\n'.format(compound_namespace))
        for compound in compounds:
            gzip_file.write('{},"{}",VENDOR;REAL\n'.
                            format(compound.c, compound.s))


def augment_original_nodes(directory, filename, has_header):
    """Augments the original nodes file and writes the relationships
    for nodes in this file to the Vendor nodes.

    :param directory: The output directory
    :param filename: The name of the original data file (to be augmented)
    :param has_header: True if the original file has a header
    """

    global num_nodes
    global num_nodes_augmented
    global num_compound_relationships

    logger.info('Augmenting {} as...'.format(filename))

    # Augmented file
    augmented_filename =\
        os.path.join(directory,
                     '{}-augmented-{}.gz'.format(output_filename_prefix,
                                                 os.path.basename(filename)))
    gzip_ai_file = gzip.open(augmented_filename, 'wt')
    # Frag to Vendor Compound relationships file
    augmented_relationships_filename =\
        os.path.join(directory,
                     '{}-molecule-compound-edges.csv.gz'.
                     format(output_filename_prefix))
    gzip_cr_file = gzip.open(augmented_relationships_filename, 'wt')
    gzip_cr_file.write(':START_ID({}),'
                       ':END_ID({}),'
                       ':TYPE\n'.format(smiles_namespace, compound_namespace))

    logger.info(' {}'.format(augmented_filename))
    logger.info(' {}'.format(augmented_relationships_filename))

    with open(filename) as i_file:

        if has_header:
            # Copy first line (header)
            hdr = i_file.readline()
            gzip_ai_file.write(hdr)

        for line in i_file:

            num_nodes += 1
            # Search for a potential MolPort identity
            # Get the MolPort compound
            # if we know the compound add a label
            augmented = False
            match_ob = enamine_re.findall(line)
            if match_ob:
                # Look for compounds where we have a costed vendor.
                # If there is one, add the label.
                for compound_id in match_ob:
                    if compound_id in compound_map:
                        new_line = line.strip() + ';V_E\n'
                        gzip_ai_file.write(new_line)
                        augmented = True
                        num_nodes_augmented += 1
                        break
                if augmented:
                    # If we've augmented the line
                    # append a relationship to the relationships file
                    # for each compound that was found...
                    for compound_id in match_ob:
                        if compound_id in compound_map:
                            # Now add vendor relationships to this row
                            frag_id = line.split(',')[0]
                            gzip_cr_file.write('"{}",{},HAS_VENDOR\n'.
                                               format(frag_id, compound_id))
                            num_compound_relationships += 1

            if not augmented:
                # No vendor for this line,
                # just write it out 'as-is'
                gzip_ai_file.write(line)

    # Close augmented nodes and the relationships
    gzip_ai_file.close()
    gzip_cr_file.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Processor (Enamine)')
    parser.add_argument('vendor_dir',
                        help='The Enamine vendor directory,'
                             ' containing the ".gz" files to be processed.')
    parser.add_argument('vendor_prefix',
                        help='The Enamine vendor file prefix,'
                             ' i.e. "June2018". Only files with this prefix'
                             ' in the vendor directory will be processed')
    parser.add_argument('nodes',
                        help='The nodes file to augment with the collected'
                             ' vendor data')
    parser.add_argument('output',
                        help='The output directory')
    parser.add_argument('--nodes-has-header',
                        help='Use if the nodes file has a header',
                        action='store_true')

    args = parser.parse_args()

    # Create the output directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if not os.path.isdir(args.output):
        error('output ({}) is not a directory'.format(args.output))

    # Process all the files...
    enamine_files = glob.glob('{}/{}*.gz'.format(args.vendor_dir, args.vendor_prefix))
    for enamine_file in enamine_files:
        extract_vendor(enamine_file)

    # Write the new nodes and relationships
    # and augment the original nodes file.
    if compound_map:
        write_compound_nodes(args.output, compound_nodes)
        augment_original_nodes(args.output, args.nodes, has_header=args.nodes_has_header)

    # Summary
    logger.info('{} compounds'.format(len(compound_map)))
    logger.info('{} nodes'.format(num_nodes))
    logger.info('{} augmented nodes'.format(num_nodes_augmented))
    logger.info('{} node compound relationships'.format(num_compound_relationships))
