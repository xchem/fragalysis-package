# Series of functions to parse input files
from frag.alysis.models import Object,Owner
from frag.utils.rdkit_utils import _parse_ligand_sdf, _get_c_of_mass, RDKitPh4, _get_waters, _get_water_coords,_parse_pdb
import math

def parse_ligands(input_file,input_type="sdf"):
    mols = _parse_ligand_sdf(input_file=input_file)
    # Now return them with their name and centre of mass
    c_of_mass_list = []
    for m in mols:
        c_of_mass_list.append(_get_c_of_mass(m))
    return c_of_mass_list


def parse_ligand_ph4s(input_file, input_type="sdf"):
    """
    Function to return a series of ligand based pharmacophores.
    :param input_file: the file to parse
    :param input_type: the type of ligands
    :return: the molecule based pharmacophores
    """
    rdkit_ph4 = RDKitPh4()
    mols = _parse_ligand_sdf(input_file=input_file)
    output_pharma_list = []
    for mol in mols:
        pharma_list = rdkit_ph4.generate_ph4_for_mol(rdmol=mol)
        output_pharma_list.append(pharma_list)
    return output_pharma_list

def parse_waters(input_pdbs, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of waters - return waters.
    :param input_pdb: the input PDB files
    :param input_mol: the input molecule (to use as a reference around which to limit)
    :return: tuple threes of coordinates of the waters
    """
    owner_list = []
    # First just get the waters from the file
    for input_pdb in input_pdbs:
        waters = _get_waters(open(input_pdb).readlines())
        water_coords = _get_water_coords(waters)
        out_l = []
        for water in water_coords:
            water = Object(water,"water")
            out_l.append(water)
        owner = Owner(out_l,input_pdb)
        owner_list.append(owner)
    return owner_list


def find_dist(mol_1_x,mol_1_y, mol_1_z, mol_2_x, mol_2_y, mol_2_z):
    """Function to find the square distance between two points in 3D
    Takes two len=3 tuples
    Returns a float"""
    return pow((mol_1_x-mol_2_x),2) + pow((mol_1_y-mol_2_y),2) + pow((mol_1_z-mol_2_z),2)


def _get_res_rmsds(input_res_list):
    """
    Helper function to get the RMSDs for a list of Residues.
    :param input_res_list: The list of RDKit molecules of residues
    :return: a list of lists of RMSDS.
    """
    # Calculate RMSD from corresponding
    num_res = len(input_res_list)
    tot_res_rmsd_list = []
    for i in range(num_res):
        this_res_rmsd_list = []
        for j in range(i,num_res):
            res_one = input_res_list[i]
            res_two = input_res_list[j]
            tot_dist = 0.0
            num_matches = 0
            for atom in res_one:
                atm1 = res_one[atom]
                if atom in res_two:
                    atm2 = res_two[atom]
                    # Find the distance
                    dist = find_dist(atm1[0], atm1[1], atm1[2],
                                 atm2[0], atm2[1], atm2[2])
                    tot_dist += dist
                    num_matches +=1
            # Find the mean square distance
            mean_dist = float(tot_dist) / float(num_matches)
            # Append the root mean square distance
            root_mean_sqr = math.sqrt(mean_dist)
            this_res_rmsd_list.append(root_mean_sqr)
        tot_res_rmsd_list.append(this_res_rmsd_list)
    return tot_res_rmsd_list

def get_res_atom_name(atom, conf):
    """
    Get the information for each atom
    :param atom: the input atom
    :param conf: the molecular conformer
    :return: the unqiue residue level name, the atom name and the position
    """
    res_info = atom.GetPDBResidueInfo()
    atom_pos =  conf.GetAtomPosition(atom.GetIdx())
    position = [atom_pos.x,atom_pos.y,atom_pos.z]
    identifiers =  [res_info.GetResidueNumber(),res_info.GetChainId(),
                            res_info.GetResidueName(), res_info.GetAltLoc()]
    unique_name = "_".join([x for x in identifiers if x != ' '])
    atom_name = res_info.GetName().strip()
    return unique_name, atom_name, position


def _get_res(input_data):
    """
    Get a list of residues with RDKit mols (RDMol,RES_NAME)
    :param input_data:
    :return:
    """
    out_dict = {}
    # Loop through the residues
    mol = _parse_pdb(input_data)
    conf = mol.GetConformer()
    atoms = mol.GetAtoms()
    for atom in atoms:
        # Get res_name
        res_name,atom_name,position = get_res_atom_name(atom,conf)
        if res_name in out_dict:
            out_dict[res_name][atom_name] = position
        else:
            out_dict[res_name] = {atom_name: position}
    return out_dict


def parse_residues(input_pdbs, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of PDB files of proteins.
    :param input_pdb: the input PDB files - with identifiers
    :return: a dict (Key Residue -> value list of molecules)
    """
    owner_list = []
    res_dict = {}
    for input_pdb in input_pdbs:
        this_res_dict = _get_res(input_pdb)
        for key in this_res_dict:
            if key in res_dict:
                res_dict[key].append(res_dict[key])
            else:
                res_dict[key] = [res_dict[key]]
    for res in res_dict:
        rmsd_coords = _get_res_rmsds(res_dict[res])
        out_l = []
        res = Object(rmsd_coords,res)
        out_l.append(res)
        owner = Owner(out_l,input_pdb)
        owner_list.append(owner)
    return owner_list