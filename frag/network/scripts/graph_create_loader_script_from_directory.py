#!/usr/bin/env python3
# coding=utf-8

"""A utility to create the Graph loader script from the files in a named
directory. This is generally employed when combining data sets,
where 2 or more graph databases are combined.

Alan Christie
July 2019
"""

import argparse
import logging
import os
import sys

from process_utils import write_load_script

# Configure basic logging
formatter = logging.Formatter('%(asctime)s %(levelname)s # %(message)s',
                              '%Y-%m-%dT%H:%M:%S')

out_hdlr = logging.StreamHandler()
out_hdlr.setFormatter(formatter)
out_hdlr.setLevel(logging.INFO)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(out_hdlr)


def error(msg):
    """Prints an error message and exists.

    :param msg: The message to print
    """
    logger.error('ERROR: {}'.format(msg))
    sys.exit(1)


parser = argparse.ArgumentParser('Graph Graph Loader Generator')
parser.add_argument('path', metavar='PATH', type=str,
                    help='The path to the files'
                         ' that are expected to be loaded.')

args = parser.parse_args()

if not os.path.isdir(args.path):
    error('path ({}) is not a directory'.format(args.output))

# Just construct a dictionary of nodes and edges based on the
# directory content. We add every file that contains the word
# 'nodes' or 'edges' and ends in '.gz'.
#
# If the file is 'nodes.csv.gz' or 'edges.csv.gz' we expect a header
# and instead of just adding the file to the list we also add the header.
# So, if we find 'nodes.csv.gz' we insert 'nodes-header.csv,nodes.csv.gz'
# into the corresponding list.
generated_files = {'nodes': [],
                   'edges': []}
possible_files = os.listdir(args.path)
for possible_file in possible_files:
    if possible_file.endswith('.gz') and \
            os.path.isfile(os.path.join(args.path, possible_file)):
        if 'edges' in possible_file:
            actual_file = possible_file
            if actual_file == 'edges.csv.gz':
                actual_file = 'edges-header.csv,edges.csv.gz'
            generated_files['edges'].append(actual_file)
        elif 'nodes' in possible_file:
            actual_file = possible_file
            if actual_file == 'nodes.csv.gz':
                actual_file = 'nodes-header.csv,nodes.csv.gz'
            generated_files['nodes'].append(actual_file)

if not generated_files['nodes']:
    error('Could not find any nodes files.')
if not generated_files['edges']:
    error('Could not find any edges files.')

logger.info('Writing load script...')
write_load_script(args.path, generated_files)
