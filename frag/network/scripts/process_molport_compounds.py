#!/usr/bin/env python

"""process_molport_compounds.py

Processes MolPort vendor compound files, expected to contain pricing
information, against the colated graph processing files.

Four new files are generated and the original nodes file
augmented with a "V_MP" label.

Note:   This module does expect `colate_all` to have been used on the original
        graph files to produce normalised supplier identities in the node file
        this module uses.

The purpose of this module is to create "Vendor" Compound and "Cost" nodes
and relationships to augment the DLS fragment database.
Every fragment line that has a MolPort identifier in the original data set
is labelled and a relationship created between it and the Vendor's compound(s).
The compounds are also related to purchasing costs for those compounds in
various "pack sizes".

Some vendor compound nodes may have no defined costs and some compounds may
not exist in the original data set.

The files generated (in a named output directory) are:

-   "molport-cost-nodes.csv.gz"
    containing nodes that define the unique set of costs.

-   "molport-compound-nodes.csv.gz"
    containing all the nodes for the vendor compounds
    (that have at least one set of pricing information).

-   "molport-compound-cost-edges.csv.gz"
    containing the "Compound" to "Cost"
    relationships using the the type of "COSTS".

-   "molport-molecule-compound-edges.csv.gz"
    containing the relationships between the original node entries and
    the "Compound" nodes. There is a relationship for every MolPort
    compound that was found in the earlier processing.

The module augments the original nodes by adding the label
"V_MP" for all MolPort compounds that have been found
to the augmented copy of the original node file that it creates.

If the original nodes file is "nodes.csv.gz" the augmented copy
(in the named output directory) will be called
"molport-augmented-nodes.csv.gz".

Note:   At the moment the original nodes.csv.gz file is expected to contain
        standardised (uppercase) compound identifiers,
        i.e. "MOLPORT:NNN-NNN-NNN" whereas the compound files (that include
        pricing information) are expected to use the supplier's identifier,
        i.e. "MolPort-NNN-NNN-NNN". See the 'molport_re' and 'supplier_prefix'
        variables in this module.

        In the future the standardiser should produce a new compound file
        that contains all the relevant columns passed through (at the moment
        it just contains SSMILES, OSMILES and ID columns.

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
cost_1mg_col = 5
cost_5mg_col = 6
cost_50mg_col = 7
blt_col = 8
expected_input_cols = {smiles_col: 'SMILES',
                       compound_col: 'MOLPORTID',
                       cost_1mg_col: 'PRICERANGE_1MG',
                       cost_5mg_col: 'PRICERANGE_5MG',
                       cost_50mg_col: 'PRICERANGE_50MG',
                       blt_col: 'BEST_LEAD_TIME'}

# The Vendor Compound node has...
# a compound id (unique for a given vendor)
# a SMILES string
# a best lead time
CompoundNode = namedtuple('CompoundNode', 'c s blt')
# The Cost node has...
# a unique id (assigned after collection)
# a pack size (mg)
# a minimum price
# a maximum price
CostNode = namedtuple('CostNode', 'ps min max')
# A unique set of cost nodes
cost_nodes = set()
# The map of Vendor Compound nodes against an array of Cost nodes
compound_cost_map = {}
# The vendor compound IDs that have pricing information
costed_compounds = set()
# All the vendor compound IDs
vendor_compounds = set()

# Prefix for output files
output_filename_prefix = 'molport'
# The namespaces of the various indices
smiles_namespace = 'F2'
compound_namespace = 'VMP'
cost_namespace = 'CMP'

# Regular expression to find the MolPort compound IDs
# (in the original nodes file).
molport_re = re.compile(r'MOLPORT:(\d+-\d+-\d+)')
supplier_prefix = 'MolPort-'

# Various diagnostic counts
num_compounds_without_costs = 0
num_compound_cost_relationships = 0
num_nodes = 0
num_nodes_augmented = 0
num_compound_relationships = 0


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


def create_cost_node(pack_size, field_value):
    """Creates a CostNode namedtuple for the provided pack size
    and corresponding pricing field. If the pricing field
    is empty or does not correspond to a recognised format
    or has no min or max value no CostNode is created.

    :param pack_size: The pack size (mg). Typically 1, 5, 50 etc.
    :param field_value: The pricing field value, e.g. "100 - 500"
    :returns: A CostNode namedtuple (or None if no pricing). The global
              set (cost_nodes) is also added to.
    """

    global cost_nodes

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
        cost_nodes.add(c_node)

    return c_node


def extract_vendor_compounds(gzip_filename):
    """Process the given file and extract vendor (and pricing) information.
    Vendor nodes are only created when there is at least one
    column of pricing information.

    :param gzip_filename: The compressed file to process
    """

    global compound_cost_map
    global costed_compounds
    global num_compounds_without_costs

    logger.info('Processing {}...'.format(gzip_filename))

    num_lines = 0
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

            smiles = fields[smiles_col]
            compound_id = fields[compound_col].split(supplier_prefix)[1]
            blt = int(fields[blt_col].strip())

            # Add the compound (a UUID) to our set of all compounds.
            # The compound ID has to be unique
            if compound_id in vendor_compounds:
                error('Duplicate compound ID ({})'.format(compound_id))
            vendor_compounds.add(compound_id)

            # Create a vendor node for this compound
            compound_node = CompoundNode(compound_id, smiles, blt)
            # Collect costs (there may be none)
            cost_node_1 = create_cost_node(1, fields[cost_1mg_col])
            cost_node_5 = create_cost_node(5, fields[cost_5mg_col])
            cost_node_50 = create_cost_node(50, fields[cost_50mg_col])

            costs = []
            if cost_node_1:
                costs.append(cost_node_1)
            if cost_node_5:
                costs.append(cost_node_5)
            if cost_node_50:
                costs.append(cost_node_50)

            # The compound cost map is a map of all compounds
            # with associated costs (which might be empty)
            compound_cost_map[compound_node] = costs

            if costs:
                # This compound has some costs,
                # add it to the costed compounds set for fast lookup later.
                costed_compounds.add(compound_id)
            else:
                num_compounds_without_costs += 1


def write_cost_nodes(directory, costs):
    """Writes the CostNodes to a node file, including a header.

    :param directory: The sub-directory to write to
    :param costs: The map of costs against their assigned UUID
    """

    filename = os.path.join(directory,
                            '{}-cost-nodes.csv.gz'.
                            format(output_filename_prefix))
    logger.info('Writing {}...'.format(filename))

    with gzip.open(filename, 'wb') as gzip_file:
        gzip_file.write(':ID({}),'
                        'currency,'
                        'pack_size_mg:INT,'
                        'min_price:FLOAT,'
                        'max_price:FLOAT,'
                        ':LABEL\n'.format(cost_namespace))
        for cost in costs:
            # Handle no value (None) in min and max,
            # replacing with an empty string...
            min = ''
            if cost.min is not None:
                min = cost.min
            max = ''
            if cost.max is not None:
                max = cost.max
            gzip_file.write('{},USD,{},{},{},COST\n'.format(costs[cost],
                                                            cost.ps,
                                                            min,
                                                            max))


def write_compound_nodes(directory, compound_cost_map):
    """Writes the CompoundNodes to a node file, including a header.

    :param directory: The sub-directory to write to
    :param compound_cost_map: The map of costs (against any available costs)
    """

    filename = os.path.join(directory,
                            '{}-compound-nodes.csv.gz'.
                            format(output_filename_prefix))
    logger.info('Writing {}...'.format(filename))

    with gzip.open(filename, 'wb') as gzip_file:
        gzip_file.write('cmpd_id:ID({}),'
                        'smiles,'
                        'best_lead_time:INT,'
                        ':LABEL\n'.format(compound_namespace))
        for compound in compound_cost_map:
            gzip_file.write('{},"{}",{},VENDOR;MOLPORT\n'.
                            format(compound.c,
                                   compound.s,
                                   compound.blt))


def write_compound_cost_relationships(directory, compound_cost_map, cost_map):
    """Writes the Vendor Compound to Costs relationships.

    :param directory: The output directory to write to
    :param compound_cost_map: A map of CompoundNodes
                              against an array of CostNodes
                              (which might be empty)
    :param cost_map: A map of CostNode to its assigned unique ID
    """

    global num_compound_cost_relationships

    filename = os.path.join(directory,
                            '{}-compound-cost-edges.csv.gz'.
                            format(output_filename_prefix))
    logger.info('Writing {}...'.format(filename))

    with gzip.open(filename, 'wb') as gzip_file:
        gzip_file.write(':START_ID({}),'
                        ':END_ID({}),'
                        ':TYPE\n'.format(compound_namespace, cost_namespace))
        for compound in compound_cost_map:
            # Generate a relationship for each cost for the compound.
            # The source is the vendor compound
            # and the destination is the Cost UUID (auto-assigned)
            for cost in compound_cost_map[compound]:
                cost_uuid = cost_map[cost]
                gzip_file.write('{},{},COSTS\n'.format(compound.c, cost_uuid))
                num_compound_cost_relationships += 1


def augment_original_nodes(directory, filename, has_header):
    """Augments the original nodes file and writes the relationships
    for nodes in this file to the Vendor nodes.
    """

    global num_nodes
    global num_nodes_augmented
    global num_compound_relationships

    logger.info('Augmenting {} as...'.format(filename))

    # Augmented file
    augmented_filename = \
        os.path.join(directory,
                     '{}-augmented-{}'.format(output_filename_prefix,
                                              os.path.basename(filename)))
    gzip_ai_file = gzip.open(augmented_filename, 'wt')
    # Frag to Vendor Compound relationships file
    augmented_relationships_filename = \
        os.path.join(directory,
                     '{}-molecule-compound-edges.csv.gz'.
                     format(output_filename_prefix))
    gzip_cr_file = gzip.open(augmented_relationships_filename, 'wt')
    gzip_cr_file.write(':START_ID({}),'
                       ':END_ID({}),'
                       ':TYPE\n'.format(smiles_namespace, compound_namespace))

    logger.info(' {}'.format(augmented_filename))
    logger.info(' {}'.format(augmented_relationships_filename))

    with gzip.open(filename, 'rt') as gzip_i_file:

        if has_header:
            # Copy first line (header)
            hdr = gzip_i_file.readline()
            gzip_ai_file.write(hdr)

        for line in gzip_i_file:

            num_nodes += 1
            # Search for a potential MolPort identity
            # Get the MolPort compound
            # if we know the compound add a label
            augmented = False
            match_ob = molport_re.findall(line)
            if match_ob:
                # Look for compounds where we have a costed vendor.
                # If there is one, add the "V_MP" label.
                for compound_id in match_ob:
                    if compound_id in costed_compounds:
                        new_line = line.strip() + ';V_MP\n'
                        gzip_ai_file.write(new_line)
                        augmented = True
                        num_nodes_augmented += 1
                        break
                if augmented:
                    # If we've augmented the line
                    # append a relationship to the relationships file
                    # for each compound that was found...
                    for compound_id in match_ob:
                        if compound_id in costed_compounds:
                            # Now add vendor relationships to this row
                            frag_id = line.split(',')[0]
                            gzip_cr_file.write('"{}",{},HAS_VENDOR\n'.format(frag_id,
                                                                             compound_id))
                            num_compound_relationships += 1

            if not augmented:
                # No vendor for this line,
                # just write it out 'as-is'
                gzip_ai_file.write(line)

    # Close augmented nodes and the relationships
    gzip_ai_file.close()
    gzip_cr_file.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Vendor Compound Processor (MolPort)')
    parser.add_argument('dir',
                        help='The MolPort vendor directory,'
                             ' containing the ".gz" files to be processed.'
                             ' All the ".gz" files in the supplied directory'
                             ' will be inspected.')
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
    molport_files = glob.glob('{}/*.gz'.format(args.dir))
    for molport_file in molport_files:
        extract_vendor_compounds(molport_file)

    # Assign unique identities to the collected Cost nodes
    # using a map of cost node against the assigned ID.
    # We have to do this now because namedtuples are immutable,
    # so we collect Cost nodes first and then create unique IDs.
    cost_uuid_map = {}
    next_cost_node_uuid = 1
    for cost_node in cost_nodes:
        cost_uuid_map[cost_node] = next_cost_node_uuid
        next_cost_node_uuid += 1

    # Write the new nodes and relationships
    # and augment the original nodes file.
    if cost_uuid_map:
        write_cost_nodes(args.output, cost_uuid_map)
    if compound_cost_map:
        write_compound_nodes(args.output, compound_cost_map)
        write_compound_cost_relationships(args.output, compound_cost_map, cost_uuid_map)
        augment_original_nodes(args.output, args.nodes, has_header=args.nodes_has_header)

    # Summary
    logger.info('{} costs'.format(len(cost_uuid_map)))
    logger.info('{} compounds with costs'.format(len(compound_cost_map)))
    logger.info('{} compounds without any costs'.format(num_compounds_without_costs))
    logger.info('{} compound cost relationships'.format(num_compound_cost_relationships))
    logger.info('{} nodes'.format(num_nodes))
    logger.info('{} augmented nodes'.format(num_nodes_augmented))
    logger.info('{} node compound relationships'.format(num_compound_relationships))
