#!/usr/bin/env python
# coding=utf-8

"""A utility to split (shatter) non-compressed (or compressed) CSV files across
'N' output CSV files where lines read from the input files with the same hash
will be placed in the same output file.

The purpose is to allow distributed deduplication as duplicate lines are
guaranteed to co-exists in the same (smaller) output files and 'sort'
can be run in a distributed fashion and the results simply concatenated.

Alan Christie
August 2019
"""

import argparse
import glob
import gzip
import os


def shatter(input_dir, input_suffix, num_files, output_basename,
            recursive=False, delete_input=False):
    """Given a list of filenames this utility places lines with the same
    hash into the same file.

    :param input_dir: The directory to find the input files
    :param input_suffix: Input file suffix (how each filename ends)
    :param num_files: The number of files to shatter to
    :param output_basename: The basename of the output file (i.e. 'nodes')
    :param recursive: True to search recursively
    :param delete_input: Deletes input files as they're scattered
    """
    assert os.path.exists(input_dir)
    assert os.path.isdir(input_dir)

    output_files = []
    for file_id in range(0, num_files):
        output_filename = '%s-%03d.csv' % (output_basename, file_id + 1)
        output_files.append(open(output_filename, 'wt'))

    # Recursive search for input files?
    # Because we're in Python 2 we can't use the simpler Python 3 glob...
    input_files = []
    if recursive:
        for root, dirs, files in os.walk(input_dir):
            for i_file in files:
                if i_file.endswith(input_suffix):
                    input_files.append(os.path.join(root, i_file))
    else:
        input_files = glob.glob('%s/*%s' % (input_dir, input_suffix))

    # For each file, scatter the lines amongst the uncompressed output files
    # based on the hash of each line (coping with .gz files if required)
    for input_file in input_files:
        if input_file.endswith('.gz'):
            with gzip.open(input_file, 'rt') as i_file:
                for line in i_file:
                    file_index = hash(line) % num_files
                    output_files[file_index].write(line)
        else:
            with open(input_file, 'rt') as i_file:
                for line in i_file:
                    file_index = hash(line) % num_files
                    output_files[file_index].write(line)
        if delete_input:
            os.remove(input_file)

    # Close all the output files
    for output_file in output_files:
        output_file.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Shatter lines amongst'
                                     ' hash-determined CSV files')
    parser.add_argument('inputDir', type=str,
                        help='The input directory')
    parser.add_argument('suffix', type=str,
                        help='The file suffix,'
                             ' Only files with this suffix'
                             ' will be processed (i.e. "nodes.csv")')
    parser.add_argument('numFiles', type=int,
                        help='The number of files to shatter to')
    parser.add_argument('outputBasename', type=str,
                        help='The basename for the output files')
    parser.add_argument("--delete-input", action='store_true',
                        help='Delete input files as they are scattered')
    parser.add_argument("--recursive", action='store_true',
                        help='Search for files recursively')

    args = parser.parse_args()
    shatter(args.inputDir, args.suffix, args.numFiles, args.outputBasename,
            args.recursive, args.delete_input)
