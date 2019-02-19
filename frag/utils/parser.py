# Series of functions to parse input files
from collections import namedtuple

from frag.utils.standardise_utils import get_standard_items, verify_header
from frag.alysis.models import Object, Owner
from frag.utils.rdkit_utils import (
    _parse_ligand_sdf,
    _get_c_of_mass,
    RDKitPh4,
    RDKitAtom,
    _get_water_coords,
    _get_waters,
    _get_res,
    _get_res_rmsds,
    _parse_pdb,
)
from rdkit import Chem

Standard = namedtuple('standard', 'smiles cmpd_id')

def _get_c_of_mass_list(mols):
    c_of_mass_list = []
    for m in mols:
        c_of_mass_list.append(_get_c_of_mass(m))
    return c_of_mass_list


def parse_ligands(input_file, input_type="sdf"):
    mols = _parse_ligand_sdf(input_file=input_file)
    # Now return them with their name and centre of mass
    return _get_c_of_mass_list(mols)


def parse_ligand_ph4s(input_mols):
    """
    Function to return a series of ligand based pharmacophores.
    :param input_mols: the RDKit molecules
    :return: the molecule based pharmacophores
    """
    rdkit_ph4 = RDKitPh4()
    rdkit_atom = RDKitAtom()
    output_pharma_list = []
    for mol in input_mols:
        if not mol:
            pharma_list = []
        else:
            pharma_list = rdkit_ph4.generate_ph4_for_mol(rdmol=mol)
            atom_list = rdkit_atom.generate_atoms_for_mol(mol)
            x, y, z = _get_c_of_mass(rdmol=mol)
            c_of_m_feat = (x, y, z, "c_of_m")
            pharma_list.append(c_of_m_feat)
            pharma_list.extend(atom_list)
        output_pharma_list.append(pharma_list)
    return output_pharma_list


def parse_waters(input_pdbs, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of waters - return waters.
    :param input_pdb: the input PDB files
    :param input_mol: the input molecule (to use as a reference around which to limit)
    :return: tuple threes of coordinates of the waters
    """
    output_water_list = []
    # First just get the waters from the file
    for input_pdb in input_pdbs:
        waters = _get_waters(open(input_pdb).readlines())
        output_water_list.append(_get_water_coords(waters))
    return output_water_list


def parse_residues(input_pdbs, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of PDB files of proteins.
    :param input_pdb: the input PDB files - with identifiers
    :return: a dict (Key Residue -> value list of molecules)
    """
    owner_list = []
    res_dict = {}
    for input_pdb in input_pdbs:
        # Loop through the residues
        mol = _parse_pdb(input_pdb)
        this_res_dict = _get_res(mol)
        for key in this_res_dict:
            if key in res_dict:
                res_dict[key].append(res_dict[key])
            else:
                res_dict[key] = [res_dict[key]]
    for res in res_dict:
        rmsd_coords = _get_res_rmsds(res_dict[res])
        out_l = []
        res = Object(rmsd_coords, res)
        out_l.append(res)
        owner = Owner(out_l, input_pdb)
        owner_list.append(owner)
    return owner_list


def get_file(file_path, output_format, file_counter):
    if output_format == "smi":
        return Chem.SmilesWriter(file_path + "_" + str(file_counter) + ".smi")
    else:
        return Chem.SDWriter(file_path + "_" + str(file_counter) + ".sdf")


def parse_mols(input_file, input_format):
    if input_format == "smi":
        return Chem.SmilesMolSupplier(input_file)
    else:
        return Chem.SDMolSupplier(input_file)


def parse_standard_file(input_file,
                        limit=0,
                        min_hac=0,
                        max_hac=0,
                        iso_flag=True):
    """Parses an Informatics Matters 'standard' SMILES file.
    The file is not expected to be compressed but is expected to contain
    columns for osmiles, isomeric and non-isomeric representations along with
    a compound identifier.

    :param input_file: The name of the standard file (expected to be compressed)
    :param limit: If non zero (+ve), limit content to no more than the
                  provided value. If used in conjunction with min/max HAC
                  then the limit will be applied to the
                  number that satisfy the HAC range will be returned,
                  rather than just the first N in the file.
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
        num_collected = 0
        for line in standard_file:

            std = get_standard_items(line)

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
