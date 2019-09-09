# coding=utf-8

from collections import namedtuple
import gzip

# The columns *every* standard file is expected to contain.
# Use UPPER_CASE.
# All standard files must start with these columns.
STANDARD_COLUMNS = ['OSMILES',
                    'ISO_SMILES',
                    'NONISO_SMILES',
                    'HAC',
                    'CMPD_ID']

# A named-tuple representation of the standard contents of
# a file line, returned by the convenient function 'get_standard_items()'.
StandardRow = namedtuple('StandardRow', 'osmiles iso noniso hac cmpd_id')
Standard = namedtuple('standard', 'smiles cmpd_id')


def verify_header(hdr_line):
    """Given the header of a standard file, this method
    raises an exception if it is not valid.

    :param hdr_line: The standard file header
    """

    line_items = hdr_line.split('\t')
    if len(line_items) < len(STANDARD_COLUMNS):
        raise Exception('Header has too few fields')
    for index in range(len(STANDARD_COLUMNS)):
        if line_items[index].strip().upper() != STANDARD_COLUMNS[index]:
            raise Exception('Expected column %s but found %s',
                            STANDARD_COLUMNS[index], line_items[index])

    # OK if we get here...


def get_standard_items(line):
    """Given a file line (that has been split), this module returns a
    named tuple of the line's standard content. Numerical values are
    converted accordingly.

    :param line: A line from a standard file.

    :returns: A named tuple or an exception if the content was in error
    """

    line_items = line.split('\t')

    # Line expected to have all our standard items.
    min_items = len(STANDARD_COLUMNS)
    num_items = len(line_items)
    if num_items < min_items:
        raise Exception('Items list is too short. Expected %d got %d',
                        min_items, num_items)

    osmiles = line_items[0].strip()
    iso = line_items[1].strip()
    noniso = line_items[2].strip()
    hac_str = line_items[3].strip()
    cmpd_id = line_items[4].strip()

    # HAC should be an integer...
    try:
        hac = int(hac_str)
    except ValueError:
        raise Exception('HAC (%s) is not an integer', hac_str)

    return StandardRow(osmiles, iso, noniso, hac, cmpd_id)


def parse_standard_file(input_file,
                        limit=0,
                        skip=0,
                        min_hac=0,
                        max_hac=0,
                        iso_flag=True):
    """Parses an Informatics Matters 'standard' SMILES file.
    The file is not expected to be compressed but is expected to contain
    columns for osmiles, isomeric and non-isomeric representations along with
    a compound identifier.

    :param input_file: The name of the standard file (expected to be uncompressed)
    :param limit: If non zero (+ve), limit content to no more than the
                  provided value. If used in conjunction with min/max HAC
                  then the limit will be applied to the
                  number that satisfy the HAC range
                  rather than just the first N molecules.
    :param skip: If non zero (+ve), skip this number of molecules before
                 considering any.
    :param min_hac: Only molecules with at least the provided number
                    of heavy atoms will be considered.
    :param max_hac: if grater than zero then only molecules with no more
                    than the provided number of heavy atoms will be considered.
    :param iso_flag: True to use the standard isomeric representation,
                     False to use the non-isomeric representation.

    :returns: a set of 'Standard' namedtuples
    """
    standards = set()
    with open(input_file, 'r') as standard_file:

        # Read (and verify) the header...
        hdr = standard_file.readline()
        verify_header(hdr)

        # Process the rest of the file...
        num_skipped = 0
        num_collected = 0
        for line in standard_file:

            std = get_standard_items(line)

            # Do we need to skip molecules before processing?
            if num_skipped < skip:
                num_skipped += 1
                continue

            # HAC within range?
            # If not, skip this line.
            if std.hac < min_hac or max_hac > 0 and std.hac > max_hac:
                continue

            # Collect..
            if iso_flag:
                standards.add(Standard(std.iso, std.cmpd_id))
            else:
                standards.add(Standard(std.noniso, std.cmpd_id))

            # Enough?
            num_collected += 1
            if limit and num_collected >= limit:
                break

    return standards


def filter_standard_file(input_file,
                         limit=0,
                         skip=0,
                         min_hac=0,
                         max_hac=0,
                         filter_file=None,
                         reject_file=None):
    """Parses an Informatics Matters 'standard' SMILES file.
    The input is filtered and written to stdout. If files are presented
    the results are written to the respective files, with the accepted
    (filtered) results going to the filter file and the rejected
    results to the reject file.

    :param input_file: The name of the standard file (expected to be compressed)
    :param limit: If non zero (+ve), limit content to no more than the
                  provided value. If used in conjunction with min/max HAC
                  then the limit will be applied to the
                  number that satisfy the HAC range
                  rather than just the first N molecules.
    :param skip: If non zero (+ve), skip this number of molecules before
                 considering any.
    :param min_hac: Only molecules with at least the provided number
                    of heavy atoms will be considered.
    :param max_hac: if grater than zero then only molecules with no more
                    than the provided number of heavy atoms will be considered.
    :param filter_file: If specified this is expected to be the name of an
                        output file to write the filtered (acceptable) results
                        to (e.g. one that ends '.tag.gz')
    :param reject_file: If specified this is expected to be the name of an
                        output file to write the rejected results to
                        (e.g. one that ends '.tag.gz'). If skip and limit are
                        also used the reject file only contains molecules
                        that were considered but rejected by the min/max HAC
                        setting. i.e. any skipped or otherwise unread lines
                        are not in the reject file.
    """
    filtered_file = None
    if filter_file:
        filtered_file = gzip.open(filter_file, 'wt')

    rejected_file = None
    if reject_file:
        rejected_file = gzip.open(reject_file, 'wt')

    with gzip.open(input_file, 'rt') as standard_file:

        # Read (and verify) the header...
        hdr = standard_file.readline()
        verify_header(hdr)
        if filtered_file:
            filtered_file.write(hdr)
        else:
            print(hdr.strip())

        if rejected_file:
            rejected_file.write(hdr)

        # Process the rest of the file...
        num_skipped = 0
        num_collected = 0
        for line in standard_file:

            std = get_standard_items(line)

            # Do we need to skip molecules before processing?
            if num_skipped < skip:
                num_skipped += 1
                continue

            # HAC within range?
            # If not, skip this line.
            if std.hac < min_hac or max_hac > 0 and std.hac > max_hac:
                # Put the 'rejected' line in a reject file?
                if rejected_file:
                    rejected_file.write(line)
                continue

            if filtered_file:
                filtered_file.write(line)
            else:
                print(line.strip())

            # Enough?
            num_collected += 1
            if limit and num_collected >= limit:
                break

    if filtered_file:
        filtered_file.close()
    if rejected_file:
        rejected_file.close()
