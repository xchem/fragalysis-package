# Tools for merging, growing and linking fragments into hybrid queries
from frag.utils.rdkit_utils import RDKitPh4, RDKitAtom


def combine_ph4(
    mols, max_num_ph4s=4, min_cluster_size=2, smooth_param=0.5, smooth_grad=0.1
):
    """
    Combine a series of pharmacophores into a single query
    :param mols: a list of RDKit molecules
    :param max_num_ph4s: the maximum number of ph4s in a query
    :param min_cluster_size: the minimum size of a cluster for a point to be allowed
    :param smooth_param: the max allowed pertubation for an atom
    :param smooth_grad: the smoothing gradient for an atom
    :return: a series of possible RDKit Pharmacophore queries
    """
    rdkit_ph4 = RDKitPh4()
    molecule_ph4s = []
    for i, mol in enumerate(mols):
        ph4s = rdkit_ph4.generate_ph4_for_mol(mol)
        molecule_ph4s.append([i, mol, ph4s])
    # Cluster and remove points that don't fit into a given cluster

    # Combine these combinatorially and generate duplicates with the smoothing


def combine_atoms(
    mols, clash_param=1.5, min_cluster_size=2, smooth_param=0.5, smooth_grad=0.1
):
    """
    Combine a series of molecules into a single list of atoms - removing clashes
    :param mols: a list of RDKit molecules
    :param clash_param: the min distance between atoms
    :param min_cluster_size: the minimum size of a cluster for a point to be allowed
    :param smooth_param: the max allowed pertubation for an atom
    :param smooth_grad: the smoothing gradient for an atom
    :return: a series of possible combinations of atoms - e.g. as input into USRCAT
    """
    rdkit_atom = RDKitAtom()
    molecule_atoms = []
    for i, mol in enumerate(mols):
        atoms = rdkit_atom.generate_atoms_for_mol(mol)
        molecule_atoms.append([i, mol, atoms])


def bond_atoms(atom_list):
    """
    Attempt to merge allowed atoms given allowed rules and combinations.
    This will be based off an R group of a given fragment.
    Loss function to define a good molecule - More bonds, more interactions, less atoms.
    Simply as ratios:
    bonds-fulfilled / atom
    interactions / atom
    Need to keep the ratio as low
    Just train on single molecules, then try merging two, then more and more.
    :param atom_list: the list of allowed atoms, x,y,z and atom desc (atomic number and hybridization state)
    :return:
    """
    pass


def generate_exclusions(proteins):
    """
    Generate a series of exclusion spheres for some cavities
    :param proteins: a list of protein cavities
    :return: the combined exclusion spheres
    """
    pass
