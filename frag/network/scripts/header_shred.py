#!/usr/bin/env python

"""
A file splitter.

A split-like utility given an input file, an output basename and the number
of lines required in each split output file. Here the input file is expected
to have a one-line header, which is replicated in each output file.

Alan Christie
February 2019
"""

import argparse
import gzip


def _split(input_stream, output_base, output_size, limit, extension):
    """Splits the lines in an input file (including a header)
    into a series of output files.

    :param input_stream: The input stream
    :type input_stream: ``stream``
    :param output_base: The basename of the output files. The files are written
                        to the current working directory with a numerical
                        suffix and .smi. If the base is `x` the first output
                        file will be named `x_1[extension]`.
    :type output_base: ``str``
    :param output_size: The number of lines in each output file.
    :type output_size: ``int``
    :param limit: Limit the total number of molecules to this value.
                  Process all if zero.
    :type limit: ``int``
    :param extension: The extension for the output files.
    :type extension: ``str``
    """

    # Get the input file's header
    header = input_stream.readline()

    file_line_count = 0
    file_number = 1
    num_molecules = 0
    output_file = None
    line = input_stream.readline()
    while line and line.strip():

        if file_line_count == 0:
            # Start a new file and write the header
            name = output_base + "_" + str(file_number) + extension
            output_file = open(name, "w")
            output_file.write(header)
            file_line_count = 0

        output_file.write(line)
        file_line_count += 1

        if file_line_count == output_size:
            # Close file
            output_file.close()
            output_file = None
            file_line_count = 0
            file_number += 1

        # Enough?
        num_molecules += 1
        if limit > 0 and num_molecules >= limit:
            break

        # Next line...
        line = input_stream.readline()

    # Done last input line.
    # Do we have anything un-written?
    if output_file:
        output_file.close()


def header_split(input_file, output_base, output_size, limit, extension=".smi"):
    """Splits the lines in an input file (including a header)
    into a series of output files.

    :param input_file: The input file
    :type input_file: ``str``
    :param output_base: The basename of the output files. The files are written
                        to the current working directory with a numerical
                        suffix and .smi. If the base is `x` the first output
                        file will be named `x_1[extension]`.
    :type output_base: ``str``
    :param output_size: The number of lines in each output file.
    :type output_size: ``int``
    :param limit: Limit the total number of molecules to this value.
                  Process all if zero.
    :type limit: ``int``
    :param extension: The extension for the output files.
    :type extension: ``str``
    """

    if input_file.endswith('.gz'):
        with gzip.open(input_file, 'rt') as smiles_file:
            _split(smiles_file, output_base, output_size, limit, extension)
    else:
        with open(input_file) as smiles_file:
            _split(smiles_file, output_base, output_size, limit, extension)


def main():
    PARSER = argparse.ArgumentParser(
        description="File splitter."
        " Given a file this utility"
        " splits the lines"
        " replicating its header"
        " to new files."
    )
    PARSER.add_argument(
        "-i", "--input_file",
        help="The name of the input file. If the input file ends '.gz'"
             " it is assumed to be compressed and a gzip reader is used.",
        required=True
    )
    PARSER.add_argument(
        "-o",
        "--output_base",
        help="The basename you want to use for your output"
        " files. Each output file will use this as a base"
        " name and will append a unique decimal number to"
        " the end before adding a .smi extension.",
        required=True,
    )
    PARSER.add_argument(
        "-s",
        "--output_size",
        help="The (maximum) number of lines in each output" " file.",
        required=True,
        type=int,
    )
    parser.add_argument('-l',
                        '--limit',
                        help='Limit processing to the first N molecules,'
                             ' process all otherwise.',
                        type=int,
                        default=0)
    ARGS = PARSER.parse_args()

    header_split(ARGS.input_file, ARGS.output_base,
                 int(ARGS.output_size), int(ARGS.limit))


if __name__ == "__main__":
    main()
