#!/usr/bin/env python
# coding=utf-8

"""standardise_utils.py

Utils shared with the various standardise modules.

The module provides a method for rendering a vendor (original) SMILES
representation into a standard form, including attributes like the
heavy atom count. As each caller may augment the standard with their own
columns then they are required to construct the output file, and its
header.

The standard (mandatory) columns are available to all in the
STANDARD_COLUMNS list. This must form the first line of every file written
and they must be the first columns in the file.

Alan Christie
February 2019
"""

import logging
from collections import namedtuple

from rdkit import Chem
from frag.utils.rdkit_utils import standardize

# The columns *every* standard file is expected to contain.
# Use UPPER_CASE.
# All standard files must start with these columns.
STANDARD_COLUMNS = ['OSMILES',
                    'ISO_SMILES',
                    'NONISO_SMILES',
                    'HAC',
                    'CMPD_ID']

# The tuple returned by calls to 'standardise()'.
# If the first field (std) is None then the standard cannot be used.
# All items are rendered as strings.
StandardInfo = namedtuple('StandardInfo', 'std iso noniso hac')
# A named-tuple representation of the standard contents of
# a file line, returned by the convenient function 'get_standard_items()'.
StandardRow = namedtuple('StandardRow', 'osmiles iso noniso hac cmpd_id')

# Out logger
logger = None

def standardise(osmiles):
    """Given a vendor (original) SMILES this method standardises
    it into a canonical form and returns a namedtuple that contains
    the standard form, the isomeric form, the non-isomeric form
    the heavy atom count (hac).

    :param osmiles: The original (non-standard) SMILES

    :return: A namedtuple containing the standard molecule
             representations and info. Errors are logged and, on error,
             the standard form will be returned as None.
    """
    global logger

    # Create our logger if it does not exist
    if not logger:
        logger = logging.getLogger(__name__)

    # Standardise and update global maps...
    # And try and handle and report any catastrophic errors
    # from dependent modules/functions.

    std = None
    iso = None
    noniso = None
    hac = 0

    mol = None
    try:
        mol = Chem.MolFromSmiles(osmiles)
    except Exception as e:
        logger.warning('MolFromSmiles(%s) exception: "%s"',
                       osmiles, e.message)
    if not mol:
        logger.error('Got nothing from MolFromSmiles(%s).'
                     ' Skipping this Vendor compound', osmiles)

    if mol:

        # Got a molecule.
        #
        # Get the HAC and try to (safely) standardise,
        # and create isomeric an non-isomeric representations.

        hac = mol.GetNumHeavyAtoms()
        try:
            std = standardize(mol)
        except Exception as e:
            logger.warning('standardize(%s) exception: "%s"',
                           osmiles, e.message)
        if not std:
            logger.error('Got nothing from standardize(%s).'
                         ' Skipping this Vendor compound', osmiles)

    if std:

        # We have a standard representation,
        # Try to generate the isomeric version...

        try:
            iso = Chem.MolToSmiles(std, isomericSmiles=True, canonical=True)
        except Exception as e:
            logger.warning('MolToSmiles(%s, iso) exception: "%s"',
                           osmiles, e.message)
        if not iso:
            logger.error('Got nothing from MolToSmiles(%s, iso).'
                         ' Skipping this Vendor compound', osmiles)

    if std:

        # We have a standard representation,
        # Try to generate the non-isomeric version...

        try:
            noniso = Chem.MolToSmiles(std, isomericSmiles=False, canonical=True)
        except Exception as e:
            logger.warning('MolToSmiles(%s, noniso) exception: "%s"',
                           osmiles, e.message)
        if not noniso:
            logger.error('Got nothing from MolToSmiles(%s, noniso).'
                         ' Skipping this Vendor compound', osmiles)

    # If anything went wrong, set std to None.
    # It's "all-or-nothing"...
    if not iso or not noniso:
        std = None

    return StandardInfo(std, iso, noniso, str(hac))


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
