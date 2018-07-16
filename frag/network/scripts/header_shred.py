#!/usr/bin/env python

"""
A file splitter.

A split-like utility given an input file, an output basename and the number
of lines required in each split output file. Here the input file is expected
to have a one-line header, which is replicated in each output file.

Alan Christie
July 2018
"""

import argparse


def header_slit(input_file, output_base, output_size, extension=".smi"):
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
    :param extension: The extension for the output files.
    :type extension: ``str``
    """

    with open(input_file) as smiles_file:

        # Get the input file's header
        smiles_file.seek(0, 0)
        header = smiles_file.readline()

        file_line_count = 0
        file_number = 1
        output_file = None
        line = smiles_file.readline()
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

            # Next line...
            line = smiles_file.readline()

        # Done last input line.
        # Do we have anything un-written?
        if output_file:
            output_file.close()


def main():
    PARSER = argparse.ArgumentParser(
        description="File splitter."
        " Given a file this utility"
        " splits the lines"
        " replicating its header"
        " to new files."
    )
    PARSER.add_argument(
        "-i", "--input_file", help="The name of the input file.", required=True
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
    ARGS = PARSER.parse_args()

    header_slit(ARGS.input_file, ARGS.output_base, int(ARGS.output_size))


if __name__ == "__main__":
    main()
