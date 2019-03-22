#!/usr/bin/env python
# coding=utf-8

"""Utilities shared by the various "process" modules.

Alan Christie
January 2019
"""

from collections import namedtuple
import gzip
import logging
import os
import sys

# Configure basic logging
logger = logging.getLogger('utils')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The Vendor Assay node has...
# assay.name
# assay.description
# assay.type
AssayNode = namedtuple('AssayNode', 'name description type')

# Augmenting report rate...
AUGMENT_REPORT_RATE = 10000000

# The 'load_neo4j.sh' script content.
# this is used in write_load_script()
# which takes a dictionary of nodes and edges
# (relationships) and formats the
# database, nodes and relationships blocks...
LOAD_SCRIPT_CONTENT = """
#!/usr/bin/env bash

# If the destination database exists
# then do nothing...
if [ ! -d /neo4j/graph/databases/{database} ]; then
    echo "Running as $(id)"
    echo "(load_neo4j.sh) $(date) Importing into '{database}'..."

    cd /data-loader
    /var/lib/neo4j/bin/neo4j-admin import \\
        --database {database} \\
        {nodes}{relationships}

    #echo "(load_neo4j.sh) Indexing..."
    #cd /var/lib/neo4j
    #/data-loader/index_neo4j.sh
    #echo "(load_neo4j.sh) Done."

    echo "(load_neo4j.sh) $(date) Imported."
fi
"""

def error(msg):
    """Prints an error message and exits.

    :param msg: The message to print
    """
    logger.error('ERROR: %s', msg)
    sys.exit(1)


def write_assay_nodes(directory,
                      output_prefix,
                      generated_files,
                      assays,
                      assay_namespace_id):
    """Writes the Assay nodes file, including a header.

    :param directory: The sub-directory to write to
    :param output_prefix: The output filename prefix
    :param generated_files: A dictionary with 'nodes' and 'edges'
                            keys that we write to when we open a file
                            for writing.
    :param assays: The set of assays to write
    :param assay_namespace_id: The indexing ID for the nodes
    """

    filename = os.path.join(directory,
                            '{}-assay-nodes.csv.gz'.
                            format(output_prefix))
    logger.info('Writing %s...', filename)

    num_nodes = 0
    generated_files['nodes'].append(filename)
    with gzip.open(filename, 'wb') as gzip_file:
        gzip_file.write('name:ID({}),'
                        'description:STRING,'
                        'type:STRING,'
                        ':LABEL\n'.format(assay_namespace_id))
        for assay in assays:
            # Write the row (omitting the trailing ';'...
            gzip_file.write('"{}","{}","{}",Activity\n'.
                            format(assay.name, assay.description, assay.type))
            num_nodes += 1

    logger.info(' {:,}'.format(num_nodes))


def write_isomol_nodes(directory,
                       output_prefix,
                       generated_files,
                       isomol_smiles,
                       isomol_namespace_id,
                       supplier_id):
    """Writes the IsoMol nodes file, including a header.

    :param directory: The sub-directory to write to
    :param output_prefix: The output filename prefix
    :param generated_files: A dictionary with 'nodes' and 'edges'
                            keys that we write to when we open a file
                            for writing.
    :param isomol_smiles: A map of standard SMILES to a list of compounds
    :param isomol_namespace_id: The graph namespace ID for the isomol record
    :param supplier_id: The supplier
    """

    filename = os.path.join(directory,
                            '{}-isomol-nodes.csv.gz'.
                            format(output_prefix))
    logger.info('Writing %s...', filename)

    num_nodes = 0
    generated_files['nodes'].append(filename)
    with gzip.open(filename, 'wb') as gzip_file:
        gzip_file.write('smiles:ID({}),'
                        'cmpd_ids:STRING[],'
                        ':LABEL\n'.format(isomol_namespace_id))
        for smiles in isomol_smiles:
            # Construct the 'array' of compounds (';'-separated)
            compound_ids_str = ''
            for compound_id in isomol_smiles[smiles]:
                compound_ids_str += '{};'.format(compound_id)
            # Write the row (omitting the trailing ';'...
            gzip_file.write('"{}",{},CanSmi;Mol;{}\n'.
                            format(smiles, compound_ids_str[:-1], supplier_id))
            num_nodes += 1

    logger.info(' {:,}'.format(num_nodes))


def write_supplier_nodes(directory,
                         output_prefix,
                         generated_files,
                         supplier_id,
                         supplier_namespace_id):
    """Writes the IsoMol nodes file, including a header.

    :param directory: The sub-directory to write to
    :param output_prefix: The output filename prefix
    :param generated_files: A dictionary with 'nodes' and 'edges'
                            keys that we write to when we open a file
                            for writing.
    :param supplier_id: The supplier
    :param supplier_namespace_id: The graph namespace ID for the supplier record
    """

    filename = os.path.join(directory,
                            '{}-supplier-nodes.csv.gz'.
                            format(output_prefix))
    logger.info('Writing %s...', filename)

    generated_files['nodes'].append(filename)
    with gzip.open(filename, 'wb') as gzip_file:
        gzip_file.write('name:ID({}),'
                        ':LABEL\n'.format(supplier_namespace_id))
        # Write the solitary row...
        gzip_file.write('"{}",Supplier\n'.format(supplier_id))

    logger.info(' 1')


def write_isomol_suppliermol_relationships(directory,
                                           output_prefix,
                                           generated_files,
                                           isomol_smiles,
                                           isomol_namespace_id,
                                           supplier_namespace_id,
                                           assay_name=None,
                                           assay_namespace_id=None,
                                           assay_compound_values=None):
    """Writes the IsoMol to SupplierMol relationships file, including a header.

    :param directory: The sub-directory to write to
    :param output_prefix: The output filename prefix
    :param generated_files: A dictionary with 'nodes' and 'edges'
                            keys that we write to when we open a file
                            for writing.
    :param isomol_smiles: The map of standardised SMILES
                          to a list of Vendor compound IDs
    :param isomol_namespace_id: The graph namespace ID for the isomol record
    :param supplier_namespace_id: The graph namespace ID for the supplier record
    """

    filename = os.path.join(directory,
                            '{}-isomol-suppliermol-edges.csv.gz'.
                            format(output_prefix))
    logger.info('Writing %s...', filename)

    gzip_iar_file = None
    if assay_name:

        # Fragment to Assay relationships file
        # Fragment to Assay relationships file
        molecule_assay_relationships_filename = \
            os.path.join(directory,
                         '{}-isomol-assay-edges.csv.gz'.
                         format(output_prefix))
        generated_files['edges'].append(molecule_assay_relationships_filename)
        gzip_iar_file = gzip.open(molecule_assay_relationships_filename, 'wt')
        gzip_iar_file.write(':START_ID({}),'
                            'value:FLOAT,'
                            ':END_ID({}),'
                            ':TYPE\n'.format(isomol_namespace_id,
                                             assay_namespace_id))

    num_edges = 0
    generated_files['edges'].append(filename)
    with gzip.open(filename, 'wb') as gzip_file:
        gzip_file.write(':START_ID({}),'
                        ':END_ID({}),'
                        ':TYPE\n'.format(isomol_namespace_id,
                                         supplier_namespace_id))
        for smiles in isomol_smiles:
            for compound_id in isomol_smiles[smiles]:

                gzip_file.write('"{}",{},HasVendor\n'.format(smiles, compound_id))
                num_edges += 1

                # IsoMol to Assay?
                if gzip_iar_file:
                    gzip_iar_file.write('"{}",{},{},Activity\n'.
                                        format(smiles,
                                               assay_compound_values[compound_id],
                                               assay_name))

    if gzip_iar_file:
        gzip_iar_file.close()

    logger.info(' {:,}'.format(num_edges))


def write_nodes(input_nodes,
                output_dir,
                output_prefix,
                generated_files,
                frag_namespace_id,
                isomol_namespace_id,
                suppliermol_namespace_id,
                isomol_smiles,
                non_isomol_isomol_smiles,
                non_isomol_smiles,
                vendor_code,
                assay_name=None,
                assay_namespace_id=None,
                assay_compound_values=None):
    """Augments the original nodes file and writes the relationships
    for nodes in this file to the Vendor nodes.

    :param generated_files: A dictionary with 'nodes' and 'edges'
                            keys that we write to when we open a file
                            for writing.
    """

    num_nodes = 0
    num_nodes_augmented = 0
    num_compound_relationships = 0
    num_compound_iso_relationships = 0

    filename = input_nodes
    logger.info('Augmenting %s as...', filename)

    # Augmented file
    augmented_filename = \
        os.path.join(output_dir,
                     '{}-augmented-{}'.format(output_prefix,
                                              os.path.basename(filename)))
    generated_files['nodes'].append(augmented_filename)
    gzip_ai_file = gzip.open(augmented_filename, 'wt')

    # Frag to SupplierMol relationships file
    augmented_noniso_relationships_filename = \
        os.path.join(output_dir,
                     '{}-molecule-suppliermol-edges.csv.gz'.
                     format(output_prefix))
    generated_files['edges'].append(augmented_noniso_relationships_filename)
    gzip_smr_file = gzip.open(augmented_noniso_relationships_filename, 'wt')
    gzip_smr_file.write(':START_ID({}),'
                        ':END_ID({}),'
                        ':TYPE\n'.format(frag_namespace_id, suppliermol_namespace_id))

    # IsoMol to Frag relationships file
    augmented_iso_relationships_filename = \
        os.path.join(output_dir,
                     '{}-isomol-molecule-edges.csv.gz'.
                     format(output_prefix))
    generated_files['edges'].append(augmented_iso_relationships_filename)
    gzip_ifr_file = gzip.open(augmented_iso_relationships_filename, 'wt')
    gzip_ifr_file.write(':START_ID({}),'
                        ':END_ID({}),'
                        ':TYPE\n'.format(isomol_namespace_id, frag_namespace_id))

    gzip_mar_file = None
    if assay_name:

        # Fragment to Assay relationships file
        molecule_assay_relationships_filename = \
            os.path.join(output_dir,
                         '{}-molecule-assay-edges.csv.gz'.
                         format(output_prefix))
        generated_files['edges'].append(molecule_assay_relationships_filename)
        gzip_mar_file = gzip.open(molecule_assay_relationships_filename, 'wt')
        gzip_mar_file.write(':START_ID({}),'
                            'value:FLOAT,'
                            ':END_ID({}),'
                            ':TYPE\n'.format(frag_namespace_id,
                                             assay_namespace_id))

    logger.info(' %s', augmented_filename)
    logger.info(' %s', augmented_noniso_relationships_filename)
    logger.info(' %s', augmented_iso_relationships_filename)

    with gzip.open(filename, 'rt') as n_file:

        # Copy (and sanity-check) the nodes header...
        expected_hdr = 'smiles:ID(%s),' \
                       'hac:INT,' \
                       'chac:INT,' \
                       'osmiles,' \
                       'cmpd_ids:STRING[],' \
                       ':LABEL' % frag_namespace_id
        hdr = n_file.readline().strip()
        if hdr != expected_hdr:
            error('Node file header is wrong, expected "{}" got "{}"'.
                  format(expected_hdr, hdr))
        gzip_ai_file.write(hdr + '\n')

        for line in n_file:

            num_nodes += 1
            # Give user a gentle reminder to stdout
            # that all is progressing...
            if num_nodes % AUGMENT_REPORT_RATE == 0:
                logger.info(' ...at node {:,} ({:,}/{:,})'.
                            format(num_nodes,
                                   num_compound_relationships,
                                   num_compound_iso_relationships))

            # Check the fragments's SMILES against our nonisomol map.
            # This is a map into our IsoMol table and is a surrogate
            # for the lack of isomeric compound IDs that the colate
            # utility should insert (but doesn't).
            #
            # Then search for MolPort compound identities on the line.

            # Get the line items
            items = line.split(',')
            # The standardised SMILES string
            frag_smiles = items[0]
            # The growing list of compound relationships for this fragment.
            # We start with the existing list of compounds from the input...
            frag_compounds = []
            existing_compounds = items[4].strip()
            if existing_compounds:
                frag_compounds = existing_compounds.split(';')
            # This fragment's labels.
            # We start with the existing list of labels from the input...
            frag_labels = []
            existing_labels = items[5].strip()
            if existing_labels:
                frag_labels = existing_labels.split(';')
            # Flag, set if there's a compound
            augment = False

            if frag_smiles in non_isomol_isomol_smiles:

                # We've found the fragment (non-iso) SMILES in the map
                # that indicates it's a non-isomeric representation
                # of an isomer. We should augment the entry.
                for noniso_isomol_smiles in non_isomol_isomol_smiles[frag_smiles]:

                    # A relationship from IsoMol to Frag.
                    gzip_ifr_file.write('"{}","{}",NonIso\n'.
                                        format(noniso_isomol_smiles,
                                               frag_smiles))

                    # A relationship (or relationships)
                    # from Frag to SupplierMol
                    if noniso_isomol_smiles in isomol_smiles:
                        for compound_id in isomol_smiles[noniso_isomol_smiles]:
                            if compound_id not in frag_compounds:
                                gzip_smr_file.write('"{}",{},HasVendor\n'.
                                                    format(frag_smiles,
                                                           compound_id))
                                frag_compounds.append(compound_id)
                    else:
                        logger.warning('Failed to find "%s" in isomol_smiles',
                                       noniso_isomol_smiles)

                    num_compound_iso_relationships += 1
                    num_compound_relationships += 1
                    augment = True

            elif frag_smiles in non_isomol_smiles:

                # A relationship (or relationships)
                # from Frag to SupplierMol
                for compound_id in non_isomol_smiles[frag_smiles]:
                    if compound_id not in frag_compounds:
                        augment = True
                        # Fragment to Supplier
                        gzip_smr_file.write('"{}",{},HasVendor\n'.
                                            format(frag_smiles,
                                                   compound_id))
                        frag_compounds.append(compound_id)
                        num_compound_relationships += 1
                        # Fragment (molecule) to Assay?
                        if gzip_mar_file:
                            gzip_mar_file.write('"{}",{},{},Activity\n'.
                                                format(frag_smiles,
                                                       assay_compound_values[compound_id],
                                                       assay_name))

            if augment:
                if 'CanSmi' not in frag_labels:
                    frag_labels.append('CanSmi')
                if 'Mol' not in frag_labels:
                    frag_labels.append('Mol')
                if vendor_code not in frag_labels:
                    frag_labels.append(vendor_code)

            # Augment the fragment line's first 4 items
            # and then adding the new compounds and label lists
            output_items = items[:4]
            output_items.append(';'.join(frag_compounds))
            output_items.append(';'.join(frag_labels))

            gzip_ai_file.write(','.join(output_items) + '\n')
            num_nodes_augmented += 1

    # Close augmented nodes and the relationships
    gzip_ai_file.close()
    gzip_smr_file.close()
    gzip_ifr_file.close()
    # And assay relations (if used)
    if gzip_mar_file:
        gzip_mar_file.close()


    return  num_nodes, \
            num_nodes_augmented, \
            num_compound_relationships, \
            num_compound_iso_relationships


def write_load_script(output_dir, generated_files):
    """Writes a neo4j-compliant load script for the files
    we generated. A dictionary of file lists indexed by
    'nodes' and 'edges' if provided.

    :param output_dir: The output directory
    :param generated_files: A dictionary of filename lists
                            keyed by 'nodes' or 'edges'
    """

    # Only expecting 'nodes' and 'edges'.
    # DO this to trap a typo where there may be
    # something like 'nodess' that would otherwise be lost
    assert 'nodes' in generated_files
    assert 'edges' in generated_files
    assert len(generated_files) == 2

    # Generate a map of variables and values to
    # that will be used to modify the script content
    script_variables = {'database': 'graph.db'}
    nodes = ''
    for entry in generated_files['nodes']:
        if nodes:
            nodes += '        '
        nodes += '--nodes "%s" \\\n' % entry
    relationships = ''
    for entry in generated_files['edges']:
        relationships += '        '
        relationships += '--relationships "%s" \\\n' % entry
    # Remove the trailing whitespace and character (\)
    relationships = relationships.rstrip()[:-1]
    script_variables['nodes'] = nodes
    script_variables['relationships'] = relationships
    # Render the file content
    load_script_content = LOAD_SCRIPT_CONTENT.format(**script_variables)
    # And write...
    script_filename = os.path.join(output_dir, 'load_neo4j.sh')
    with open(script_filename, 'w') as script_file:
        script_file.write(load_script_content.strip())
