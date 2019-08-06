#!/usr/bin/env python3
# coding=utf-8

"""A utility to split (shatter) non-compressed CSV files in N-CSV files where
lines with the same hash will be placed in the same file. The purpose
is to allow distributed deduplication as duplicate lines are guaranteed to
co-exists in the same (smaller) files.

Alan Christie
August 2019
"""

import argparse
import glob
import os


def shatter(input_dir, input_suffix, num_files, output_basename):
    """Given a list of filenames this utility places lines with the same
    hash into the same file.

    :param input_dir: The directory to find the input files
    :param input_suffix: Input file suffix (how each filename ends)
    :param num_files: The number of files to shatter to
    :param output_basename: The basename of the output file (i.e. 'nodes')
    """

    output_files = []
    for file_id in range(0, num_files):
        output_filename = '%s-%03d.csv' % (output_basename, file_id + 1)
        output_files.append(open(output_filename, 'wt'))

    assert os.path.exists(input_dir)
    assert os.path.isdir(input_dir)
    input_files = glob.glob('%s/*%s' % (input_dir, input_suffix))
    for input_file in input_files:
        with open(input_file, 'rt') as i_file:
            for line in i_file:
                file_index = hash(line) % num_files
                output_files[file_index].write(line)

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

    args = parser.parse_args()
    shatter(args.inputDir, args.suffix, args.numFiles, args.outputBasename)
