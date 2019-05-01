#!/usr/bin/env python3
# coding=utf-8

"""A utility to create the Graph loader script from the files in a named
directory.

Alan Christie
May 2019
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
# 'nodes' or 'edges'.
generated_files = {'nodes': [],
                   'edges': []}
possible_files = os.listdir(args.path)
for possible_file in possible_files:
    if 'edges' in possible_file:
        file_path = os.path.join(args.path, possible_file)
        if os.path.isfile(file_path):
            generated_files['edges'].append(possible_file)
    elif 'nodes' in possible_file:
        file_path = os.path.join(args.path, possible_file)
        if os.path.isfile(file_path):
            generated_files['nodes'].append(possible_file)

if not generated_files['nodes']:
    error('Could not find any nodes files.')
if not generated_files['edges']:
    error('Could not find any edges files.')

logger.info('Writing load script...')
write_load_script(args.path, generated_files)
