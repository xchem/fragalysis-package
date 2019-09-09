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

# The tuple returned by calls to 'standardise()'.
# If the first field (std) is None then the standard cannot be used.
# All items are rendered as strings.
StandardInfo = namedtuple('StandardInfo', 'std iso noniso hac')

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
