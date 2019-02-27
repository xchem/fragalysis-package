#!/usr/bin/env python
# coding=utf-8

"""process_molport_compounds.py

Processes standardised MolPort vendor compound files,
expected to contain pricing information,
against the prepared (CSV) graph files.

For the graph design refer to the Google-Drive graph model document at...

    https://drive.google.com/file/d/1g4jT3yhwQYqsKwMpE3fYAA7dgGBYhBiw

The files generated (in a named output directory) are:

-   "molport-suppliermol-nodes.csv.gz"
    containing nodes that define the unique set of compound IDs
    and their original smiles.

-   "molport-supplier-nodes.csv.gz"
    containing the supplier node(s) for the vendor.

-   "molport-suppliermol-supplier-edges.csv.gz"
    containing the "SupplierMol" to "Supplier"
    relationships using the the type of "Availability".

Every fragment line that has a MolPort identifier in the original data set
is labelled and a relationship created between it and the Vendor's compound(s).
The compounds are also related to purchasing costs for those compounds in
various "pack sizes".

Some vendor compound nodes may have no defined costs and some compounds may
not exist in the original data set.

The files generated (in a named output directory) are:

-   "molport-isomol-nodes.csv.gz"
    containing information about compounds that are isomeric.

-   "molport-molecule-suppliermol-edges.csv.gz"
    containing the relationships between the original fragment node entries
    and the Vendor "Compound" nodes (where the compounds not isomeric)

-   "molport-isomol-molecule-edges.csv.gz"
    containing the relationships between IsoMol entries and
    the fragment nodes (where the fragment is isomeric).

The module augments the original nodes by adding the located compound IDs,
and a labels for all MolPort compounds that have been
found in the original node file.

If the original nodes file is "nodes.csv.gz" the augmented copy
(in the named output directory) will be called
"molport-augmented-nodes.csv.gz".

-   "molport-unknown-fragment-compounds.txt"
    is a file that contains vendor compounds referred to in the fragment file
    that are not in the Vendor data.

Alan Christie
January 2019
"""

import argparse
from collections import namedtuple
import gzip
import logging
import os
import re
import sys

from process_utils import error
from process_utils import write_supplier_nodes
from process_utils import write_isomol_nodes
from process_utils import write_isomol_suppliermol_relationships
from process_utils import write_nodes

# Configure basic logging
logger = logging.getLogger('molport')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The minimum number of columns in the input data (a standardised file).
# Essentially a map of expected column names indexed by column number.
expected_min_num_cols = 8
osmiles_col = 0
iso_smiles_col = 1
noniso_smiles_col = 2
hac_col = 3
compound_col = 4
cost_col = {1: 4, 5: 5, 50: 6}
blt_col = 7
expected_input_cols = {osmiles_col: 'OSMILES',
                       iso_smiles_col: 'ISO_SMILES',
                       noniso_smiles_col: 'NONISO_SMILES',
                       hac_col: 'HAC',
                       compound_col: 'CMPD_ID',
                       cost_col[1]: 'PRICERANGE_1MG',
                       cost_col[5]: 'PRICERANGE_5MG',
                       cost_col[50]: 'PRICERANGE_50MG',
                       blt_col: 'BEST_LEAD_TIME'}

# The Vendor SupplierMol node has...
# a compound id (unique for a given vendor)
# a SMILES string
SupplierMolNode = namedtuple('SupplierMol', 'cmpd_id osmiles')
# The Vendor Supplier node has...
# a name
SupplierNode = namedtuple('Supplier', 'name')
# The Cost node has...
# a pack size (mg)
# a minimum price
# a maximum price
CostNode = namedtuple('CostNode', 'ps min max')

# Map of Vendor compounds that are isomeric, and their standard representation.
# The index is a Vendor compound ID and the value is the standardised form.
# If the compound is in this map it is isometric.
compound_isomer_map = {}
# Map of standardised SMILES to vendor compound(s)
# that have isomeric representations.
# The index is standardised (isomeric) SMILES
# and the value is a set() of Vendor compound IDs
isomol_smiles = {}
# Map of standardised SMILES to vendor compounds(s)
# that are not isomeric representations.
non_isomol_smiles = {}
# Map of non-isomeric SMILES representations to isomeric smiles
# (where the molecule is isomeric). This helps lookup
# Vendor molecules that are isomeric rather than using the
# Vendor's compound ID.
non_isomol_isomol_smiles = {}
# All the vendor compound IDs
vendor_compounds = set()
# The set of all vendor compounds found in the fragment line
# where a Vendor compound was not found.
unknown_vendor_compounds = set()

# The supplier symbolic name
supplier_name = 'MolPort'
# Prefix for output files
output_filename_prefix = 'molport'
# The namespaces of the various indices
frag_namespace = 'F2'
suppliermol_namespace = 'SM_MP'
supplier_namespace = 'S'
isomol_namespace = 'ISO'

# Regular expression to find the MolPort compound IDs
# (in the original nodes file).
molport_re = re.compile(r'MOLPORT:(\d+-\d+-\d+)')
# The compound identifier prefix
# the vendor uses in the the compound files...
supplier_prefix = 'MolPort-'
# The prefix we use in our fragment file
# and the prefix we use for our copy of the
molport_prefix = 'MOLPORT:'

# Various diagnostic counts
num_vendor_iso_mols = 0
num_vendor_mols = 0
num_vendor_molecule_failures = 0

# The line rate at which the augmenter writes updates to stdout.
# Every 20 million?
augment_report_rate = 20000000


def create_cost_node(pack_size, field_value):
    """Creates a CostNode namedtuple for the provided pack size
    and corresponding pricing field. If the pricing field
    is empty or does not correspond to a recognised format
    or has no min or max value no CostNode is created.

    :param pack_size: The pack size (mg). Typically 1, 5, 50 etc.
    :param field_value: The pricing field value, e.g. "100 - 500"
    :returns: A CostNode namedtuple (or None if no pricing)
    """

    # The cost/pricing field value
    # has a value that is one of:
    #
    # "min - max"   e.g. "50 - 100"
    # "< max"       e.g. "< 1000"
    # "> min"       e.g. "> 50"

    min_val = None
    max_val = None
    c_node = None
    if field_value.startswith('>'):
        min_val = float(field_value.split()[1])
    elif field_value.startswith('<'):
        max_val = float(field_value.split()[1])
    elif ' - ' in field_value:
        min_val = float(field_value.split(' - ')[0])
        max_val = float(field_value.split(' - ')[1])

    if min_val is not None or max_val is not None:
        c_node = CostNode(pack_size, min_val, max_val)

    return c_node


def extract_vendor_compounds(suppliermol_gzip_file,
                             suppliermol_edges_gzip_file,
                             supplier_id,
                             gzip_filename,
                             limit):
    """Process the given file and extract vendor (and pricing) information.
    Vendor nodes are only created when there is at least one
    column of pricing information.

    This method extracts vendor information and writes the following files: -

    -   "molport-suppliermol-nodes.csv.gz"
    -   "molport-suppliermol-supplier-edges.csv.gz"

    The following files are expected to be written elsewhere: -

    -   "molport-supplier-nodes.csv.gz"

    The "ID" in the SupplierMol nodes file is the Compound ID and the
    "ID" of the (single) Supplier node is the supplier Name.

    As we load the Vendor compounds we 'standardise' the SMILES and
    determine whether they represent an isomer or not.

    :param suppliermol_gzip_file: The SupplierMol node file
    :param suppliermol_edges_gzip_file: The SupplierMol to Supplier edges file
    :param supplier_id: The ID of the supplier node
    :param gzip_filename: The compressed standard file to process
    :param limit: If non-zero, limit precessing to only the first N molecules

    :returns: The number of molecules processed
    """

    global compound_isomer_map
    global isomol_smiles
    global non_isomol_smiles
    global non_isomol_isomol_smiles
    global num_vendor_iso_mols
    global num_vendor_mols
    global num_vendor_molecule_failures

    logger.info('Processing %s...', gzip_filename)

    num_lines = 0
    num_processed = 0
    with gzip.open(gzip_filename, 'rt') as gzip_file:

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
            if field_names[col_num].strip() != expected_input_cols[col_num]:
                error('expected "{}" in column {} found "{}"'.
                      format(expected_input_cols[col_num],
                             col_num,
                             field_names[col_num]))

        # Columns look right...

        for line in gzip_file:

            num_lines += 1
            fields = line.split('\t')
            if len(fields) <= 1:
                continue

            osmiles = fields[osmiles_col]
            iso = fields[iso_smiles_col]
            noniso = fields[noniso_smiles_col]
            compound_id = fields[compound_col]
            blt = int(fields[blt_col].strip())

            # Add the compound (expected to be unique)
            # to our set of 'all compounds'.
            if compound_id in vendor_compounds:
                error('Duplicate compound ID ({})'.format(compound_id))
            vendor_compounds.add(compound_id)

            # Is it isomeric?
            num_vendor_mols += 1
            if iso != noniso:

                # Yes
                num_vendor_iso_mols += 1
                if iso not in isomol_smiles:
                    # This standardised SMILES is not
                    # in the map of existing isomers
                    # so start a new list of customer compounds...
                    new_set = set()
                    new_set.add(compound_id)
                    isomol_smiles[iso] = new_set
                else:
                    # Standard SMILES already
                    isomol_smiles[iso].add(compound_id)
                compound_isomer_map[compound_id] = iso
                # Put a lookup of iso representation from the non-iso
                if noniso not in non_isomol_isomol_smiles:
                    new_set = set()
                    new_set.add(iso)
                    non_isomol_isomol_smiles[noniso] = new_set
                else:
                    non_isomol_isomol_smiles[noniso].add(iso)

            else:

                # Not an isomeric representation
                if noniso not in non_isomol_smiles:
                    new_set = set()
                    new_set.add(compound_id)
                    non_isomol_smiles[noniso] = new_set
                else:
                    non_isomol_smiles[noniso].add(compound_id)

            # Write the SupplierMol entry
            suppliermol_gzip_file.write('{},"{}",Available\n'.
                                        format(compound_id,
                                               osmiles))

            # And add suitable 'Availability' relationships with the Supplier
            for quantity in [1, 5, 50]:
                cost = create_cost_node(quantity, fields[cost_col[quantity]])
                if cost:
                    cost_min = str(cost.min) if cost.min else ''
                    cost_max = str(cost.max) if cost.max else ''
                    suppliermol_edges_gzip_file.\
                        write('{},{},{},{},USD,{},{},Availability\n'.
                              format(compound_id,
                                     quantity,
                                     cost_min,
                                     cost_max,
                                     blt,
                                     supplier_id))

            # Enough?
            num_processed += 1
            if limit and num_processed >= limit:
                break

    return num_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Processor (MolPort)')
    parser.add_argument('vendor_file',
                        help='The vendor standardised file (gzipped).')
    parser.add_argument('input_nodes',
                        help='The compressed (gzipped) nodes file to'
                             ' augment with the collected vendor data')
    parser.add_argument('output',
                        help='The output directory')
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise.')

    args = parser.parse_args()

    # Create the output directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if not os.path.isdir(args.output):
        error('output ({}) is not a directory'.format(args.output))

    # -------
    # Stage 1 - Process our standardised Vendor File
    # -------

    # Open new files for writing.
    #
    # The output files are: -
    # - One for the SupplierMol nodes
    # - And one for the relationships to the (expected) supplier node.
    # - And one for the imomeric molecules.
    suppliermol_filename = os.path.\
        join(args.output,
             '{}-suppliermol-nodes.csv.gz'.
             format(output_filename_prefix))
    logger.info('Writing %s...', suppliermol_filename)
    suppliermol_gzip_file = gzip.open(suppliermol_filename, 'wt')
    suppliermol_gzip_file.write('cmpd_id:ID({}),'
                                'osmiles,'
                                ':LABEL\n'.format(suppliermol_namespace))

    suppliermol_edges_filename = os.path.\
        join(args.output,
             '{}-suppliermol-supplier-edges.csv.gz'.
             format(output_filename_prefix))
    logger.info('Writing %s...', suppliermol_edges_filename)
    suppliermol_edges_gzip_file = gzip.open(suppliermol_edges_filename, 'wt')
    suppliermol_edges_gzip_file.write(':START_ID({}),'
                                      'quantity,'
                                      'price_min,'
                                      'price_max,'
                                      'currency,'
                                      'lead_time,'
                                      ':END_ID({}),'
                                      ':TYPE\n'.format(suppliermol_namespace,
                                                       supplier_namespace))

    _ = extract_vendor_compounds(suppliermol_gzip_file,
                                 suppliermol_edges_gzip_file,
                                 supplier_name,
                                 args.vendor_file,
                                 args.limit)

    # Close the SupplierMol and the edges file.
    suppliermol_gzip_file.close()
    suppliermol_edges_gzip_file.close()

    # Write the supplier node file...
    write_supplier_nodes(args.output,
                         output_filename_prefix,
                         supplier_name,
                         supplier_namespace)

    # -------
    # Stage 2 - Write the IsoMol nodes
    # -------
    # We have collected and written SupplierMol nodes, Supplier nodes
    # and relationships and have a map of the vendor molecules
    # that are isomeric.

    write_isomol_nodes(args.output,
                       output_filename_prefix,
                       isomol_smiles,
                       isomol_namespace,
                       supplier_name)
    write_isomol_suppliermol_relationships(args.output,
                                           output_filename_prefix,
                                           isomol_smiles,
                                           isomol_namespace,
                                           supplier_namespace)

    # -------
    # Stage 3 - Augment
    # -------
    # Augment the processed nodes file
    # and attach relationships between it, the IsoMol and SupplierMol nodes.
    # This stage: -
    # - Creates up to 2 new relationships between:
    #   - IsoMol and Fragment Network
    #   - Fragment Network and SupplierMol

    num_nodes, \
    num_nodes_augmented, \
    num_compound_relationships, \
    num_compound_iso_relationships = write_nodes(args.input_nodes,
                                                 args.output,
                                                 output_filename_prefix,
                                                 frag_namespace,
                                                 isomol_namespace,
                                                 supplier_namespace,
                                                 isomol_smiles,
                                                 non_isomol_isomol_smiles,
                                                 non_isomol_smiles,
                                                 'V_MP')

    # Now complete we write a "done" file to the output.
    # Processing may be time-consuming
    # so this file helps us avoid unnecessary re-processing on failure.
    # This can be used by the automation (ansible) framework to
    # decide whether processing was completed successfully.
    # If there's a 'done' file we can safely assume that processing
    # is complete.
    open(os.path.join(args.output, 'done'), 'a').close()

    # Summary
    logger.info('{:,}/{:,} vendor molecules/iso'.format(num_vendor_mols, num_vendor_iso_mols))
    logger.info('{:,} vendor molecule failures'.format(num_vendor_molecule_failures))
    logger.info('{:,}/{:,} nodes/augmented'.format(num_nodes, num_nodes_augmented))
    logger.info('{:,}/{:,} node compound relationships/iso'.format(num_compound_relationships, num_compound_iso_relationships))

    # Dump compounds that were referenced in the fragment file
    # but not found in the vendor data.
    # Or remove any file that might already exist.
    unknown_vendor_compounds_file_name = os.path.join(args.output,
                             '{}-unknown_vendor_compounds.txt'.
                             format(output_filename_prefix))
    if unknown_vendor_compounds:
        file_name = os.path.join(args.output,
                                 '{}-unknown_vendor_compounds.txt'.
                                 format(output_filename_prefix))
        logger.info('{:,} unknown compounds (see {})'.
                    format(len(unknown_vendor_compounds),
                           unknown_vendor_compounds_file_name))
        with open(unknown_vendor_compounds_file_name, 'wt') as unknown_vendor_compounds_file:
            for unknown_vendor_compound in unknown_vendor_compounds:
                unknown_vendor_compounds_file.write(unknown_vendor_compound + '\n')
    else:
        logger.info('0 unknown compounds')
        if os.path.exists(unknown_vendor_compounds_file_name):
            os.remove(unknown_vendor_compounds_file_name)
