#!/usr/bin/env python

"""process_emolecules_plus_compounds.py

Processes standardised eMolecules vendor files, not expected to contain
pricing information.

For the graph design refer to the Google-Drive graph model document at...

    https://drive.google.com/file/d/1g4jT3yhwQYqsKwMpE3fYAA7dgGBYhBiw

The purpose of this module is to create "Vendor" Compound nodes
and relationships to augment the DLS fragment database.
Every fragment line that has an eMolecules identifier in the original data set
is labelled and a relationship created between it and the Vendor's compound(s).

Some vendor compound nodes may not exist in the original data set.

The files generated (in a named output directory) are:

-   "emolecules-plus-compound-nodes.csv.gz"
    containing all the nodes for the vendor compounds.

-   "emolecules-plus-molecule-compound_edges.csv.gz"
    containing the relationships between the original node entries and
    the "Vendor" nodes. There is a relationship for every Enamine
    compound that was found in the earlier processing.

The module augments the original nodes by adding the label
"V_EMOLS_PLUS" for all compounds that have been found
to the augmented copy of the original node file that it creates.

If the original nodes file is "nodes.csv" the augmented copy
(in the named output directory) will be called
"emolecules-plus-augmented-nodes.csv.gz".

Important Note  If the graph content changes in any way the
--------------  our graph_version must be changed.
                So, if the format of the node or relationship files change
                (e.g.  new columns, new labels, new or modified anything)
                our version must change.

Alan Christie
May 2019
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
from process_utils import write_load_script

# Configure basic logging
logger = logging.getLogger('emols')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# Our graph version.
#
# Altered if any change occurs that changes the topology of the graph.
# Its format is the date of the change ('YYYY-MM-DD') followed by
# a dot-delimited number that's incremented for each change on that day.
# i.e. '2019-05-21.2' is the second version on the 21st May 2019.
graph_version = '2019-05-26.1'

# The minimum number of columns in the input data (a standardised file).
# Essentially a map of expected column names indexed by column number.
expected_min_num_cols = 5
osmiles_col = 0
iso_smiles_col = 1
noniso_smiles_col = 2
hac_col = 3
compound_col = 4
expected_input_cols = {osmiles_col: 'OSMILES',
                       iso_smiles_col: 'ISO_SMILES',
                       noniso_smiles_col: 'NONISO_SMILES',
                       hac_col: 'HAC',
                       compound_col: 'CMPD_ID'}

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
supplier_name = 'EMOLS_PLUS'
# Prefix for output files
output_filename_prefix = 'emolecules-plus'
# The namespaces of the various indices
suppliermol_namespace = 'SM_EMOLS_PLUS'
supplier_namespace = 'S'
isomol_namespace = 'ISO-EMOLS_PLUS'
vendor_code = 'V_EMOLS_PLUS'

# The list of files generated.
# Used to generate the accompanying `load_neo4j.sh`.
# We add to this every time we open a file for writing.
#
# There is an implicit 'edges.csv.gz' and we add a header
EDGES_HDR_FILENAME = 'edges-header.csv'
generated_files = {'nodes': [],
                   'edges': ['{},edges.csv.gz'.format(EDGES_HDR_FILENAME)]}

# Various diagnostic counts
num_nodes = 0
num_nodes_augmented = 0
num_compound_relationships = 0
num_compound_iso_relationships = 0
num_vendor_iso_mols = 0
num_vendor_mols = 0
num_vendor_molecule_failures = 0


def extract_vendor_compounds(suppliermol_gzip_file,
                             suppliermol_edges_gzip_file,
                             supplier_id,
                             gzip_filename,
                             limit,
                             min_hac,
                             max_hac):
    """Process the given file and extract vendor (and pricing) information.
    Vendor nodes are only created when there is at least one
    column of pricing information.

    This method extracts vendor information and writes the following files: -

    -   "enamine-suppliermol-nodes.csv.gz"
    -   "enamine-suppliermol-supplier-edges.csv.gz"

    The following files are expected to be written elsewhere: -

    -   "enamine-supplier-nodes.csv.gz"

    The "ID" in the SupplierMol nodes file is the Compound ID and the
    "ID" of the (single) Supplier node is the supplier Name.

    As we load the Vendor compounds we 'standardise' the SMILES and
    determine whether they represent an isomer or not.

    :param suppliermol_gzip_file: The SupplierMol node file
    :param suppliermol_edges_gzip_file: The SupplierMol to Supplier edges file
    :param supplier_id: The ID of the supplier node
    :param gzip_filename: The compressed standard file to process
    :param limit: If non-zero, limit precessing to only the first N molecules
    :param min_hac: Minimum HAC (0 for no minimum)
    :param max_hac: Maximum HAC (0 for no maximum)

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
            hac = int(fields[hac_col])
            iso = fields[iso_smiles_col]
            noniso = fields[noniso_smiles_col]
            compound_id = fields[compound_col]

            # If min/max HAC have been provided
            # use them to eliminate compounds.
            if hac < min_hac:
                continue
            elif max_hac and hac > max_hac:
                continue

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

            # And add a suitable 'Availability' relationship with the Supplier
            suppliermol_edges_gzip_file. \
                write('{},{},Availability\n'.
                      format(compound_id,
                             supplier_id))

            # Enough?
            num_processed += 1
            if limit and num_processed >= limit:
                break

    return num_processed


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Processor (eMolecules)')
    parser.add_argument('vendor_file',
                        help='The vendor standardised file (gzipped).')
    parser.add_argument('input_nodes',
                        help='The compressed (gzipped) nodes file to'
                             ' augment with the collected vendor data')
    parser.add_argument('output',
                        help='The output directory')
    parser.add_argument('-r', '--replace-input',
                        dest="replace_input", action="store_true",
                        help='When processing is complete replace the'
                             ' input nodes file with the augmented output.'
                             ' If used the load script is not generated.')
    parser.add_argument('--processing-version',
                        type=str, default='undefined',
                        help='The graph processing version (fragalysis version).'
                             ' Used as a property in the Supplier node.')
    parser.add_argument('--process-id',
                        type=str, default='undefined',
                        help='The process ID (the origin of the data).'
                             ' Used as a property in the Supplier node.')
    parser.add_argument('--build-number',
                        type=int, default=0,
                        help='The anticipated build number (0 if undefined).'
                             ' Used as a property in the Supplier node.')
    parser.add_argument('-l', '--limit',
                        type=int, default=0,
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise.')
    parser.add_argument('--min-hac',
                        type=int, default=0,
                        help='Limit processing to molecules with at least this'
                             ' number of heavy atoms')
    parser.add_argument('--max-hac',
                        type=int, default=0,
                        help='Limit processing to molecules with no more than'
                             ' this number of heavy atoms')

    args = parser.parse_args()

    # Create the output directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if not os.path.isdir(args.output):
        error('output ({}) is not a directory'.format(args.output))

    # Sanity-check key arguments
    if args.build_number < 1:
        error('build-number cannot be less then 1 ({})'.format(args.build_number))
    if args.limit < 0:
        error('limit cannot be -ve ({})'.format(args.limit))
    if args.min_hac < 0:
        error('min-hac cannot be -ve ({})'.format(args.min_hac))
    if args.max_hac < 0:
        error('max-hac cannot be -ve ({})'.format(args.max_hac))

    # -------
    # Stage 1 - Process our standardised Vendor File
    # -------

    # Open new files for writing.
    #
    # The output files are: -
    # - One for the SupplierMol nodes
    # - And one for the relationships to the (expected) supplier node.
    # - And one for the imomeric molecules.
    suppliermol_filename = os.path. \
        join(args.output,
             '{}-suppliermol-nodes.csv.gz'.
             format(output_filename_prefix))
    logger.info('Writing %s...', suppliermol_filename)
    generated_files['nodes'].append(suppliermol_filename)
    suppliermol_gzip_file = gzip.open(suppliermol_filename, 'wt')
    suppliermol_gzip_file.write('cmpd_id:ID({}),'
                                'osmiles,'
                                ':LABEL\n'.format(suppliermol_namespace))

    suppliermol_edges_filename = os.path. \
        join(args.output,
             '{}-suppliermol-supplier-edges.csv.gz'.
             format(output_filename_prefix))
    logger.info('Writing %s...', suppliermol_edges_filename)
    generated_files['edges'].append(suppliermol_edges_filename)
    suppliermol_edges_gzip_file = gzip.open(suppliermol_edges_filename, 'wt')
    suppliermol_edges_gzip_file.write(':START_ID({}),'
                                      ':END_ID({}),'
                                      ':TYPE\n'.format(suppliermol_namespace,
                                                       supplier_namespace))

    _ = extract_vendor_compounds(suppliermol_gzip_file,
                                 suppliermol_edges_gzip_file,
                                 'eMolecules-BB',
                                 args.vendor_file,
                                 args.limit,
                                 args.min_hac,
                                 args.max_hac)

    # Close the SupplierMol and the edges file.
    suppliermol_gzip_file.close()
    suppliermol_edges_gzip_file.close()

    # Write the supplier node file...
    write_supplier_nodes(args.output,
                         output_filename_prefix,
                         generated_files,
                         supplier_name,
                         supplier_namespace,
                         vendor_code,
                         graph_version,
                         args.processing_version,
                         args.process_id,
                         args.build_number,
                         args.limit,
                         args.min_hac,
                         args.max_hac)

    # -------
    # Stage 2 - Write the IsoMol nodes
    # -------
    # We have collected and written SupplierMol nodes, Supplier nodes
    # and relationships and have a map of the vendor molecules
    # that are isomeric.

    write_isomol_nodes(args.output,
                       output_filename_prefix,
                       generated_files,
                       isomol_smiles,
                       isomol_namespace,
                       supplier_name)
    write_isomol_suppliermol_relationships(args.output,
                                           output_filename_prefix,
                                           generated_files,
                                           isomol_smiles,
                                           isomol_namespace,
                                           suppliermol_namespace)

    # -------
    # Stage 3 - Augment
    # -------
    # Augment the processed nodes file
    # and attach relationships between it, the IsoMol and SupplierMol nodes.
    # This stage: -
    # - Creates up to 2 new relationships between:
    #   - IsoMol and Fragment Network
    #   - Fragment Network and SupplierMol

    augmented_file_path, \
    num_nodes, \
    num_nodes_augmented, \
    num_compound_relationships, \
    num_compound_iso_relationships = write_nodes(args.input_nodes,
                                                 args.output,
                                                 output_filename_prefix,
                                                 generated_files,
                                                 isomol_namespace,
                                                 suppliermol_namespace,
                                                 isomol_smiles,
                                                 non_isomol_isomol_smiles,
                                                 non_isomol_smiles,
                                                 vendor_code)

    # Replace input file (normally used as part of a chain,
    # i.e. during the combination playbook).
    # Otherwise create the loader script here.
    # If we replace the output file the loader script is something
    # that has to be generated separately.
    if args.replace_input:
        logger.info('Replacing input (%s -> %s)...',
                    augmented_file_path, args.input_nodes)
        # Remove the original input file and repl;ace it with the output
        os.remove(args.input_nodes)
        os.rename(augmented_file_path, args.input_nodes)
    else:
        # Before we finish,
        # write a convenient loader script
        # for all the files we generated...
        logger.info('Writing load script...')
        write_load_script(args.output, generated_files)

    # Finish by writing the expected edges header file...
    edges_header_file = open(EDGES_HDR_FILENAME, 'wt')
    edges_header_file.write(':START_ID(F2),:END_ID(F2),label,:TYPE\n')
    edges_header_file.close()

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
