#!/usr/bin/env python

"""process_senp7_compounds.py

Processes standardised senp7 vendor files.

The files generated (in a named output directory) are:

-   "senp7-compound-nodes.csv.gz"
    containing all the nodes for the vendor compounds.

-   "senp7-molecule-compound_edges.csv.gz"
    containing the relationships between the original node entries and
    the "Vendor" nodes. There is a relationship for every senp7
    compound that was found in the earlier processing.

The module augments the original nodes by adding the label
"V_HTS" for all MolPort compounds that have been found
to the augmented copy of the original node file that it creates.

If the original nodes file is "nodes.csv" the augmented copy
(in the named output directory) will be called
"hts-augmented-nodes.csv.gz".

Alan Christie
November 2018
"""

import argparse
import gzip
import logging
import os
import sys

from process_utils import error
from process_utils import write_supplier_nodes
from process_utils import write_isomol_nodes
from process_utils import write_isomol_suppliermol_relationships
from process_utils import write_nodes
from process_utils import write_assay_nodes
from process_utils import AssayNode

# Configure basic logging
logger = logging.getLogger('hts')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The minimum number of columns in the input files and
# and a map of expected column names indexed by (0-based) column number.
expected_min_num_cols = 6
osmiles_col = 0
iso_smiles_col = 1
noniso_smiles_col = 2
hac_col = 3
compound_col = 4
activity_col = 5
expected_input_cols = {osmiles_col: 'OSMILES',
                       iso_smiles_col: 'ISO_SMILES',
                       noniso_smiles_col: 'NONISO_SMILES',
                       hac_col: 'HAC',
                       compound_col: 'CMPD_ID',
                       activity_col: 'INHIB_5UM'}

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
# Compound ID to Assay activity map.
assay_compound_values = {}

# The supplier symbolic name
supplier_name = 'SENP7'
# Prefix for output files
output_filename_prefix = 'senp7'
# The namespaces of the various indices
frag_namespace = 'F2'
suppliermol_namespace = 'SM_SENP7'
supplier_namespace = 'S'
isomol_namespace = 'ISO-SENP7'
assay_namespace = 'A_SENP7'

# Various diagnostic counts
num_vendor_iso_mols = 0
num_vendor_mols = 0
num_vendor_molecule_failures = 0
num_activity_value_failures = 0


def extract_vendor_compounds(suppliermol_gzip_file,
                             suppliermol_edges_gzip_file,
                             supplier_id,
                             gzip_filename,
                             limit):
    """Process the given file and extract vendor information.

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
    global assay_compound_values
    global num_vendor_iso_mols
    global num_vendor_mols
    global num_vendor_molecule_failures
    global num_activity_value_failures

    logger.info('Processing {}...'.format(gzip_filename))

    num_lines = 0
    num_processed = 0
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

            osmiles = fields[osmiles_col]
            iso = fields[iso_smiles_col]
            noniso = fields[noniso_smiles_col]
            compound_id = fields[compound_col]
            activity = 0.0
            try:
                activity = float(fields[activity_col])
            except ValueError:
                num_activity_value_failures += 1

            # The compound ID must be unique
            if compound_id in vendor_compounds:
                error('Duplicate compound "{}"'.format(compound_id))
            vendor_compounds.add(compound_id)

            # Add activity to compound-id map...
            assay_compound_values[compound_id] = activity

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

            # And add a suitable 'Availability' relationship with the Supplier
            suppliermol_edges_gzip_file. \
                write('{},{},Availability\n'.
                      format(compound_id, supplier_id))

            # Enough?
            num_processed += 1
            if limit and num_processed >= limit:
                break

    return num_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Processor (HTS/SENP7)')
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
    # - And one for the relationships to the (expected) supplier node
    # - And one for the molecule to assay relationships
    suppliermol_filename = os.path. \
        join(args.output,
             '{}-suppliermol-nodes.csv.gz'.
             format(output_filename_prefix))
    logger.info('Writing %s...', suppliermol_filename)
    suppliermol_gzip_file = gzip.open(suppliermol_filename, 'wt')
    suppliermol_gzip_file.write('cmpd_id:ID({}),'
                                'osmiles,'
                                ':LABEL\n'.format(suppliermol_namespace))

    suppliermol_edges_filename = os.path. \
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

    mol_assay_edge_filename = os.path. \
        join(args.output,
             '{}-molecule-assay-edges.csv.gz'.
             format(output_filename_prefix))
    logger.info('Writing %s...', mol_assay_edge_filename)
    mol_assay_edge_gzip_file = gzip.open(mol_assay_edge_filename, 'wt')
    mol_assay_edge_gzip_file.write(':START_ID({}),'
                                   'osmiles,'
                                   'value,'
                                   ':END_ID({}),'
                                   ':TYPE\n'.format(frag_namespace,
                                                    assay_namespace))

    _ = extract_vendor_compounds(suppliermol_gzip_file,
                                 suppliermol_edges_gzip_file,
                                 supplier_name,
                                 args.vendor_file,
                                 args.limit)

    # Close the opened files.
    suppliermol_gzip_file.close()
    suppliermol_edges_gzip_file.close()
    mol_assay_edge_gzip_file.close()

    # Write the supplier node file...
    write_supplier_nodes(args.output,
                         output_filename_prefix,
                         supplier_name,
                         supplier_namespace)

    # Write assay node file...
    # There's just one assay here.
    # Our 'assay' record...
    assay_name = 'Percent inhibition'
    assay_description = 'Percent inhibition against SENP7'
    assay_type = '%INH'

    assays = [AssayNode(assay_name, assay_description, assay_type)]
    write_assay_nodes(args.output,
                      output_filename_prefix,
                      assays,
                      assay_namespace)

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
                                           supplier_namespace,
                                           assay_name,
                                           assay_namespace,
                                           assay_compound_values)

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
                                                 'V_HTS',
                                                 assay_name,
                                                 assay_namespace,
                                                 assay_compound_values)

    # Summary
    logger.info('{:,}/{:,} vendor molecules/iso'.format(num_vendor_mols, num_vendor_iso_mols))
    logger.info('{:,} vendor molecule failures'.format(num_vendor_molecule_failures))
    logger.info('{:,} vendor activity value failures'.format(num_activity_value_failures))
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
