#!/usr/bin/env python
# coding=utf-8

"""process_prep.py

Pre-process raw nodes.txt and edges.txt files
(the output of initial graph processing)
so that they can be used in the vendor processing sequence.

This module simply renders the raw text file with a header and comma-separated
fields. The output nodes file is also pre-populated with compound ID and label
fields in the nodes file.

Alan Christie
May 2019
"""

import argparse
import gzip
import logging
import os
import sys

# Configure basic logging
logger = logging.getLogger('process-prep')
out_hdlr = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)
logger.addHandler(out_hdlr)
logger.setLevel(logging.INFO)

# The fragment node namespace
# A label in the nodes csv file.
frag_namespace = 'F2'

# The line rate at which the prop writes updates to stdout.
# Every 20 million?
prep_report_rate = 20000000


def error(msg):
    """Prints an error message and exits.

    :param msg: The message to print
    """
    logger.error('ERROR: %s', msg)
    sys.exit(1)


def prep(input_dir, output_dir):
    """Preps nodes and edges files.

    :param input_dir: The input directory (expected to contain nodes.txt
                      and edges.txt.gz)
    :param output_dir: The output directory (known to exist)
    :returns: False on failure
    """

    nodes_txt_filename = os.path.join(input_dir, 'nodes.txt.gz')
    if not os.path.isfile(nodes_txt_filename):
        error("Can't find nodes file ({})".format(nodes_txt_filename))
    edges_txt_filename = os.path.join(input_dir, 'edges.txt.gz')
    if not os.path.isfile(edges_txt_filename):
        error("Can't find edges file ({})".format(edges_txt_filename))

    # Convert the nodes file...
    # This is white-space delimited file with five fields: -
    #
    #   NODE B1OCCO1.CBr 7 5 C1CCCC1.CBr

    logger.info('Processing %s...', nodes_txt_filename)

    csv_filename = os.path.join(output_dir, 'nodes.csv.gz')
    csv_file = gzip.open(csv_filename, 'wt')
    # Start the output with a neo-suitable header...
    columns = ["smiles:ID(%s)" % frag_namespace,
               "hac:INT",
               "chac:INT",
               "osmiles",
               "cmpd_ids:STRING[]",
               ":LABEL"]
    csv_file.write(','.join(columns) + '\n')

    line_num = 1
    expected_field_count = 5
    with gzip.open(nodes_txt_filename, 'rt') as txt_file:
        for line in txt_file:
            items = line.split()

            # Sanity check the line...
            #
            # This is important because un-trapped format errors
            # here can lead to obscure problems downstream or
            # when loading data into the graph database.
            if len(items) != expected_field_count:
                logger.error('Node line %s does not have %s fields": %s',
                             line_num,
                             expected_field_count,
                             line.strip())
                return False
            # The first field must be the value 'NODE'
            if items[0] != 'NODE':
                logger.error('Node line %s does not start with "NODE": %s',
                             line_num,
                             line.strip())
                return False
            # There can only be one 'NODE' string in the line
            if line.count('NODE') > 1:
                logger.error('Node line %s has more than one "NODE" string: %s',
                             line_num,
                             line.strip())
                return False

            # Insert (empty) compound IDs and the initial label.
            # These are augmented in vendor-specific processing modules.
            items.append('')
            items.append(frag_namespace)
            # Write the items out (omitting the first column 'NODE')
            csv_file.write(','.join(items[1:]) + '\n')

            # Give user a gentle reminder to stdout
            # that all is progressing...
            if line_num % prep_report_rate == 0:
                logger.info(' ...at node {:,}'.format(line_num))
            line_num += 1

    csv_file.close()

    # Convert the edges file...
    # This is white-space delimited file with four fields: -
    #
    #   EDGE Br.C Br FG|C|C|FG|Br|Br

    logger.info('Processing %s...', edges_txt_filename)

    csv_filename = os.path.join(output_dir, 'edges.csv.gz')
    csv_file = gzip.open(csv_filename, 'wt')
    columns = [":START_ID(%s)" % frag_namespace,
               ":END_ID(%s)" % frag_namespace,
               "label",
               ":TYPE"]
    csv_file.write(','.join(columns) + '\n')

    line_num = 1
    expected_field_count = 4
    expected_label_delimiter_count = 5
    label_delimiter = '|'
    with gzip.open(edges_txt_filename, 'rt') as txt_file:
        for line in txt_file:
            items = line.split()

            # Sanity check the line...
            if len(items) != expected_field_count:
                logger.error('Edge line %s does not have %s fields": %s',
                             line_num,
                             expected_field_count,
                             line.strip())
                return False
            if items[0] != 'EDGE':
                logger.error('Edge line %s does not start with "EDGE": %s',
                             line_num,
                             line.strip())
                return False
            # There can only be one 'EDGE' string in the line
            if line.count('EDGE') > 1:
                logger.error('Edge line %s has more than one "EDGE" string: %s',
                             line_num,
                             line.strip())
                return False
            # There have to be the right number of edge label delimiters ("|")
            edge_label = items[3]
            label_delimiter_count = edge_label.count(label_delimiter)
            if label_delimiter_count != expected_label_delimiter_count:
                logger.error('Edge line %s label has wrong delimiter count: %s',
                             line_num,
                             line.strip())
                return False

            # Insert our type column value (required)
            items.append('FRAG')
            # Write the items out (omitting the first column 'EDGE')
            csv_file.write(','.join(items[1:]) + '\n')

            # Give user a gentle reminder to stdout
            # that all is progressing...
            line_num += 1
            if line_num % prep_report_rate == 0:
                logger.info(' ...at edge {:,}'.format(line_num))

    csv_file.close()

    # OK if we get here
    return True

if __name__ == '__main__':

    parser = argparse.\
        ArgumentParser('Vendor Compound Pre-Processor',
                       description='Pre-processes the nodes'
                                   ' and edges text files ready for'
                                   ' subsequent vendor augmentation')
    parser.add_argument('input',
                        help='The directory containing compressed (gzipped)'
                             ' raw text nodes and edges files (expected to'
                             ' be called nodes.txt.gz and edges.txt.gz).')
    parser.add_argument('output',
                        help='The output directory where the node and edge'
                             ' compressed csv files will be written')

    args = parser.parse_args()

    # Create the output directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if not os.path.isdir(args.output):
        error('output ({}) is not a directory'.format(args.output))

    if not os.path.isdir(args.input):
        error('Input ({}) is not a directory'.format(args.input))

    # Process but exit on error
    if not prep(args.input, args.output):
        sys.exit(1)

    # Now complete we write a "done" file to the output.
    # Processing may be time-consuming
    # so this file helps us avoid unnecessary re-processing on failure.
    # This can be used by the automation (ansible) framework to
    # decide whether processing was completed successfully.
    # If there's a 'done' file we can safely assume that processing
    # is complete.
    open(os.path.join(args.output, 'done'), 'a').close()
