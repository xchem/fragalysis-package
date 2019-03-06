#!/usr/bin/env python

"""
A file joiner.

A join-like utility given input files, an output basename and the number
of lines required in each split output file. Here the input file is expected
to have a one-line header, which is replicated in each output file.

Alan Christie
March 2019
"""

import argparse
import gzip
import os


def header_join(input_dir, prefix, output_filename):
    """Joins the lines of a series of input files (including a header)
    into a single output file with only one header.

    :param input_dir: The input directory
    :type input_dir: ``str``
    :param prefix: The input file prefix
    :type prefix: ``str``
    :param output_filename: The output file.
    :type output_filename: ``str``
    """

    # List all the files in the inout directory.
    # We'll only process those that begin with the given prefix.
    possible_files = os.listdir(input_dir)

    output_file = gzip.open(output_filename + '.gz', 'wt')
    written_header = False

    for possible_file in possible_files:
        if possible_file.startswith(prefix):
            file_path = os.path.join(input_dir, possible_file)
            if possible_file.endswith('.gz'):

                with gzip.open(file_path, 'rb') as input_stream:

                    header = input_stream.readline()
                    if not written_header:
                        output_file.write(header)
                        written_header = True

                    line = input_stream.readline()
                    while line and line.strip():
                        output_file.write(line)
                        line = input_stream.readline()

            else:

                with open(file_path) as input_stream:

                    header = input_stream.readline()
                    if not written_header:
                        output_file.write(header)
                        written_header = True

                    line = input_stream.readline()
                    while line and line.strip():
                        output_file.write(line)
                        line = input_stream.readline()

    output_file.close()

def main():
    PARSER = argparse.ArgumentParser(
        description="File joiner."
        " Given a directory this utility"
        " joins the lines replicating one header to new (gzipped) file."
    )
    PARSER.add_argument(
        "input_dir",
        help="The name of the input directory.",
    )
    PARSER.add_argument(
         "prefix",
         help="The file prefix,"
              " i.e. 'HTS_'. Only files with this prefix"
              " will be processed")
    PARSER.add_argument(
        "output_file",
        help="The filename you want to use for your output file."
             " The filename will be given the suffix '.gz'.",
    )
    ARGS = PARSER.parse_args()

    header_join(ARGS.input_dir, ARGS.prefix, ARGS.output_file)


if __name__ == "__main__":
    main()
